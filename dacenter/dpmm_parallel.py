import argparse
import copy
import random
import time
import pickle
from multiprocessing import Pool
from collections import defaultdict
import numpy as np
from scipy.special import gammaln
import matplotlib.pyplot as plt
import seaborn as sns

from dacenter_run import load_vmatrix

plt.style.use('ggplot')


class Cluster:
    def __init__(self, alpha, L):
        self.alpha = alpha   # concentration hyperparameter
        self.log_gamma = gammaln(self.alpha) - 2. * gammaln(0.5 * self.alpha)   # constant
        self.L = L

        self.base_count = np.full((2, self.L), 0., dtype=float)   # base count *without* prior
        self.ndata = 0   # cluster size
        self.log_factorial = 0.   # for calculation of Ewens

    def add(self, x):
        self.base_count[x, range(self.L)] += 1
        self.ndata += 1
        if self.ndata > 2:
            self.log_factorial += np.log(self.ndata - 1)

    def remove(self, x):
        self.base_count[x, range(self.L)] -= 1
        self.ndata -= 1
        if self.ndata > 2:
            self.log_factorial -= np.log(self.ndata)

    # Calculate posterior probability of this cluster
    def calculate_cluster_parameter(self):
        self.cluster_parameter = (self.base_count + 0.5 * self.alpha) / (np.sum(self.base_count, axis=0) + self.alpha)

    # <precomputed> == True indicates that cluster parameter is already computed
    def logp_x(self, x, precomputed=False):
        if precomputed:
            return np.sum(np.log(self.cluster_parameter[x, range(self.L)]))
        else:
            return np.sum(np.log((self.base_count[x, range(self.L)] + 0.5 * self.alpha) / (np.sum(self.base_count, axis=0) + self.alpha)))

    def calculate_likelihood(self, precomputed=False):
        logp = self.L * self.log_gamma   # TODO: is this formula OK?

        # XXX: prior should be pevious posterior?

        if not precomputed:
            self.calculate_cluster_parameter()
        logp += np.sum(self.base_count * np.log(self.cluster_parameter))
        self.likelihood = logp


class Clustering:
    def __init__(self,
                 data,
                 alpha,
                 n_parallel,
                 start_from_one_cluster=True,
                 verbose=False):
        self.x = data   # observed units
        self.alpha = alpha   # concentration hyperparameter
        self.n_parallel = n_parallel   # degree of parallel proposals in each iteration
        self.verbose = verbose

        self.N, self.L = self.x.shape   # number of units and variant sites
        self.s = np.zeros(self.N, dtype=int)   # cluster assignment of x
        self.c = {}   # clusters
        self._init_clusters(start_from_one_cluster)

        # Pairs of clusters whose merge proposal was rejected before
        # This is used to avoid redundant calculations of the acceptance rate
        # To add a pair:
        #    "rejected_merge_pairs[<cluster_index_1>].add(<cluster_index_2>)"
        #    "rejected_merge_pairs[<cluster_index_2>].add(<cluster_index_1>)"
        # When cluster(s) is updated:
        #    "del rejected_merge_pairs[<cluster_index>]"
        # To check whether a pair is registered:
        #    "if <cluster_index_1> in rejected_merge_pairs[<cluster_index_2>] and <cluster_index_2> in rejected_merge_pairs[<cluster_index_1>]"
        self.rejected_merge_pairs = defaultdict(set)

        # constant used for posterior calculation
        self.log_AF = np.sum(np.log(np.array(range(self.N)) + self.alpha))
        # constant for new cluster generation in Gibbs sampling
        self.log_normalizer_const = gammaln(self.alpha) - gammaln(self.alpha + 1) + gammaln(0.5 * self.alpha + 1) - gammaln(0.5 * self.alpha)
        #self.log_normalizer_const = np.log(0.5)   # simplified version when alpha = 1

        # initial cluster likelihood
        # initial cluster parameters are also calculated
        for cluster in list(self.c.values()):
            cluster.calculate_likelihood()

        # store posterior probabilities and # of clusters for each iteration
        # clustering with the highest posterior so far is also stored
        self.max_post = self.calculate_posterior()   # initial posterior
        self.max_s = self.s.copy()   # initial clustering
        self.posterior_ncluster = [(self.max_post, len(self.c))]

    def _init_clusters(self, start_from_one_cluster):
        if start_from_one_cluster:
            # All data belong to single cluster
            cluster = Cluster(self.alpha, self.L)
            c_idx = 0
            self.c[c_idx] = cluster
            for i in range(self.N):
                self.c[c_idx].add(self.x[i])
                self.s[i] = c_idx
        else:
            # Each data belongs to distinct cluster
            for i in range(self.N):
                cluster = Cluster(self.alpha, self.L)
                self.c[i] = cluster
                self.c[i].add(self.x[i])
                self.s[i] = i

    def show_clusters(self, c):
        for i, cluster in c.items():
            print(i, cluster.base_count)

    def calculate_posterior(self):
        # Ewens and cluster likelihood
        # NOTE: all cluster likelihoods must be updated during samplings
        logp = len(self.c) * np.log(self.alpha) - self.log_AF
        logp += sum([cluster.log_factorial + cluster.likelihood
                     for cluster in list(self.c.values())])
        return logp

    def plot_transition(self, start=0, end=None, figsize=(10, 8)):
        fig, ax1 = plt.subplots(figsize=figsize)
        ax1.plot([x[0] for x in self.posterior_ncluster[start:end]], color='r')
        ax2 = ax1.twinx()
        ax2.plot([x[1] for x in self.posterior_ncluster[start:end]], color='b')
        plt.grid(False)
        plt.show()

    # Make a new cluster with a single data
    # Cluster index is same as the data index
    def _make_new_cluster(self, idx, s, c):
        cluster = Cluster(self.alpha, self.L)
        c[idx] = cluster
        c[idx].add(self.x[idx])
        s[idx] = idx

    def _split_merge_sampling(self, arg):
        c_idx_old, init_i, init_j, x_restricted = arg
        propose_type = "split" if len(c_idx_old) == 1 else "merge"

        if propose_type == "split":   # Split
            # Calculate s and c only for new clusters
            s_launch = {}   # cluster assignment of x
            c_launch = {}   # clusters

            # Create new clusters consisting of only x[init_i] and x[init_j], respectively
            self._make_new_cluster(init_i, s_launch, c_launch)
            self._make_new_cluster(init_j, s_launch, c_launch)

            #c_restricted = [init_i, init_j]   # restricted (new) cluster indices for init_i and init_j
            c_idx_new = [init_i, init_j]

            if self.verbose:
                print("**SPLIT** cluster %d by data %d and %d" % (c_idx_old[0], init_i, init_j))

            # TODO: initial assignment by the column with highest entropy
            # that is, argmin_{col} sum(col) - len(col)/2
            # or, Monte Carlo according to the relative entropy

            ## Perform restricted Gibbs sampling (first sequential scan + several random sampling) on x_resticted
            # First sequential scan (batch update at last)
            p_cluster = np.zeros((2, len(x_restricted)))   # [:, 0] for init_i, [:, 1] for init_j
            for i, c_idx in enumerate(c_idx_new):
                c_launch[c_idx].calculate_cluster_parameter()   # precompute parameters of the new clusters (consisting of only 1 data here)
                p_cluster[i, :] = np.exp(np.array([c_launch[c_idx].logp_x(self.x[idx], precomputed=True) for idx in x_restricted]))   # single-vs-single
            p_acc_cluster = np.add.accumulate(p_cluster / np.sum(p_cluster, axis=0))

            # First assignments to one of the two clusters
            for j, idx in enumerate(x_restricted):
                rnd = random.random()
                for k, p in enumerate(p_acc_cluster[:, j]):
                    if p >= rnd:
                        c_idx = c_idx_new[k]
                        c_launch[c_idx].add(self.x[idx])
                        s_launch[idx] = c_idx
                        break

            #if verbose:
            #    print("#FIRST ASSIGNMENT BY SEQUENTIAL SCAN#")
            #    print([x[1] for x in sorted(s_launch.items(), key=lambda x: x[0])])

            # Additional Gibbs sampling
            #n_iterations = self.c[c_idx_old[0]].ndata   # same as the size of the original cluster
            n_iterations = 0   # only sequential scan
            for ii in range(n_iterations):   # TODO: improve this criterion
                # Drop one (restricted) data randomly
                r = x_restricted[random.randint(0, len(x_restricted) - 1)]
                c_launch[s_launch[r]].remove(self.x[r])

                # Calculate probabilities of assignment
                logp_cluster = np.array([c_launch[c_idx].logp_x(self.x[r]) for c_idx in c_idx_new])   # single-vs-single
                p_cluster = np.exp(logp_cluster) * np.array([c_launch[c_idx].ndata for c_idx in c_idx_new])
                p_acc_cluster = np.add.accumulate(p_cluster / np.sum(p_cluster))

                # Assignment
                rnd = random.random()
                for k, p in enumerate(p_acc_cluster):
                    if p >= rnd:
                        c_idx = c_idx_new[k]
                        c_launch[c_idx].add(self.x[r])
                        s_launch[r] = c_idx
                        break

                # TODO: should calculate posterior probability and select a result with highest value?

            if self.verbose:
                print("#LAST ASSIGNMENT AFTER GIBBS SAMPLING#")
                #print([x[1] for x in sorted(s_launch.items(), key=lambda x: x[0])])

            # This precomputed cluster parameters will be used for the calculations of both q and L below
            for c_idx in c_idx_new:
                c_launch[c_idx].calculate_cluster_parameter()   # Update cluster parameters of the new clusters here

            # Calculate q(c^split | c)
            logq = 0.
            for j in x_restricted:   # init_i and init_j are omitted   # TODO: is this proper?
                k = 0 if s_launch[j] == c_idx_new[0] else 1   # assigned cluster
                c_launch[c_idx_new[k]].remove(self.x[j])   # temporarily remove the data
                # use precomputed parameter for the cluster to which the focal data below does not belong
                logp_cluster = np.array([c_launch[c_idx].logp_x(self.x[j]) if i == k else c_launch[c_idx].logp_x(self.x[j], precomputed=True) for i, c_idx in enumerate(c_idx_new)])
                p_cluster = np.exp(logp_cluster) * np.array([c_launch[c_idx].ndata for c_idx in c_idx_new])
                logq += np.log(p_cluster[k]) - np.log(np.sum(p_cluster))
                c_launch[c_idx_new[k]].add(self.x[j])   # restore the data

            # Calculate probabilities used for acceptance ratio
            logp_transition = -logq

            logp_ewens = np.log(self.alpha)
            logp_ratio = 0.

            # new splitted clusters
            for c_idx in c_idx_new:
                logp_ewens += c_launch[c_idx].log_factorial
                cluster = c_launch[c_idx]
                cluster.calculate_likelihood(precomputed=True)   # Update likelihoods of the new clusters here
                logp_ratio += cluster.likelihood
                if self.verbose:
                    print("logp_cluster_new", c_idx, "=", cluster.likelihood)

            # old single cluster
            logp_ewens -= self.c[c_idx_old[0]].log_factorial
            cluster = self.c[c_idx_old[0]]   # likelihood of the original clusters should be already calculated
            logp_ratio -= cluster.likelihood
            if self.verbose:
                print("logp_cluster original =", cluster.likelihood)

        else:   # Merge
            if self.verbose:
                print("**MERGE** cluster %d and %d" % (c_idx_old[0], c_idx_old[1]))

            if c_idx_old[0] in self.rejected_merge_pairs[c_idx_old[1]] and c_idx_old[1] in self.rejected_merge_pairs[c_idx_old[0]]:
                # This merge proposal is previously rejected
                if self.verbose:
                    print("previously rejected")
                return (False, None, c_idx_old)

            # Set initial launch state
            c_launch = {}
            # We have to copy the original clusters because we perform sampling in them
            for c_idx in c_idx_old:
                c_launch[c_idx] = copy.deepcopy(self.c[c_idx])

            # Calculate q(c | c^merge)
            # NOTE: parameters of the clusters to be merged should be already precomputed in the posterior update
            logq = 0.
            for j in x_restricted:
                k = 0 if self.s[j] == c_idx_old[0] else 1   # originally assigned cluster
                c_launch[c_idx_old[k]].remove(self.x[j])   # temporarily remove the data
                # use precomputed parameter for the cluster to which the focal data below does not belong
                logp_cluster = np.array([c_launch[c_idx].logp_x(self.x[j]) if i == k else c_launch[c_idx].logp_x(self.x[j], precomputed=True) for i, c_idx in enumerate(c_idx_old)])
                p_cluster = np.exp(logp_cluster) * np.array([c_launch[c_idx].ndata for c_idx in c_idx_old])
                logq += np.log(p_cluster[k]) - np.log(np.sum(p_cluster))
                c_launch[c_idx_old[k]].add(self.x[j])   # restore the data

            # Merge two clusters (index c_restricted[0] will be merged to index c_restricted[1])
            c_idx_new = [c_idx_old[1]]   # list for integrality with c_idx_old
            c_launch[c_idx_new[0]].add(self.x[init_i])   # because init_i is not in x_restricted
            for j in x_restricted:
                if self.s[j] == c_idx_old[0]:
                    c_launch[c_idx_new[0]].add(self.x[j])

            # Calculate probabilities used for acceptance ratio
            logp_transition = logq

            logp_ewens = -np.log(self.alpha)
            logp_ewens += c_launch[c_idx_new[0]].log_factorial
            for j in range(2):
                logp_ewens -= self.c[c_idx_old[j]].log_factorial

            logp_ratio = 0.
            cluster = c_launch[c_idx_new[0]]
            cluster.calculate_likelihood()   # Update likelihood of the new cluster here; cluster parameter is also calculated
            if self.verbose:
                print("logp_cluster_new", c_idx_new[0], "=", cluster.likelihood)
            logp_ratio += cluster.likelihood
            for c_idx in c_idx_old:
                cluster = self.c[c_idx]
                if self.verbose:
                    print("logp_cluster original", c_idx, "=", cluster.likelihood)
                logp_ratio -= cluster.likelihood

        # Proposal
        logp_proposal = min(0, logp_transition + logp_ewens + logp_ratio)
        if self.verbose:
            print(logp_proposal, "=", logp_transition, "+", logp_ewens, "+", logp_ratio)
        p_proposal = np.exp(logp_proposal)

        rnd = random.random()
        if p_proposal >= rnd:
            accepted = True
            cluster_new = [c_launch[c_idx] for c_idx in c_idx_new]
            if propose_type == "split":
                x_split = [[idx for idx, _c_idx in s_launch.items() if _c_idx == c_idx] for c_idx in c_idx_new]
            else:
                x_split = [x_restricted + [init_i, init_j]]
            changes = (cluster_new, x_split)
            if self.verbose:
                print("accepted")
        else:
            accepted = False
            changes = None
            if self.verbose:
                print("rejected")

        return (accepted, changes, c_idx_old)

    def split_merge_sampling(self):
        # Sample two distinct data
        init_i, init_j = random.sample(range(self.N), k=2)

        # NOTE: distinguish between split and merge by "len(c_idx_old)"
        if self.s[init_i] == self.s[init_j]:   # split
            c_idx_old = [self.s[init_i]]
        else:   # merge
            c_idx_old = [self.s[init_i], self.s[init_j]]

        x_restricted = [i for i in range(self.N) if i != init_i and i != init_j and self.s[i] in c_idx_old]

        args = (c_idx_old, init_i, init_j, x_restricted)
        accepted, changes, c_idx_old = self._split_merge_sampling(args)

        max_post_updated = False
        if accepted:
            # if split:
            #    cluster_new = [new cluster 0 induced by init_i, new cluster 1 induced by init_j]
            #    x_split = [[indices of x in new cluster 0], [indices of x in new cluster 1]]
            # if merge:
            #    cluster_new = [new merged cluster]
            #    x_split = [x_restricted]
            cluster_new, x_split = changes

            # Update master <c_launch> and <s_launch>
            for i, cluster in enumerate(cluster_new):   # i = 0, 1 if split; i = 0 if merge
                c_idx = max(list(self.c.keys())) + 1
                self.c[c_idx] = cluster
                for idx in x_split[i]:   # if merge, s_changes == [x_restricted]
                    self.s[idx] = c_idx
            for c_idx in c_idx_old:
                del self.c[c_idx]

            # delete previously-rejected information
            for c_idx in c_idx_old:
                if c_idx in self.rejected_merge_pairs:
                    del self.rejected_merge_pairs[c_idx]

            # Update posterior etc if necessary
            logp_post = self.calculate_posterior()
            self.posterior_ncluster.append((logp_post, len(self.c)))
            if self.max_post < logp_post:   # max posterior is updated
                self.max_post = logp_post
                self.max_s = self.s.copy()
                max_post_updated = True
                if self.verbose:
                    print("posterior: %f -> %f" % (self.max_post, logp_post))
        else:   # rejected
            if len(c_idx_old) == 2:   # add to the rejected-pair list
                self.rejected_merge_pairs[c_idx_old[0]].add(c_idx_old[1])
                self.rejected_merge_pairs[c_idx_old[1]].add(c_idx_old[0])

            self.posterior_ncluster.append(self.posterior_ncluster[-1])

        return max_post_updated

    def split_merge_sampling_parallel(self):
        # Calculate multiple independent proposals in parallel (***OBSOLETE because multiprocessing for each iteration is too time-consuming***)
        args = []
        x_remaining = set(range(self.N))
        while len(args) < self.n_parallel and len(x_remaining) >= 2:
            # Sample two distinct data
            init_i, init_j = random.sample(x_remaining, k=2)
            x_remaining.remove(init_i)
            x_remaining.remove(init_j)

            # NOTE: distinguish between split and merge by "len(c_idx_old)"
            if self.s[init_i] == self.s[init_j]:   # split
                c_idx_old = [self.s[init_i]]
            else:   # merge
                c_idx_old = [self.s[init_i], self.s[init_j]]

            x_restricted = [i for i in list(x_remaining) if self.s[i] in c_idx_old]
            for i in x_restricted:
                x_remaining.remove(i)

            args.append((c_idx_old, init_i, init_j, x_restricted))

        # Reflect every accepted proposals on this launch state
        s_launch = self.s.copy()
        c_launch = copy.deepcopy(self.c)

        updated = False
        exe_pool = Pool(len(args))
        for accepted, changes, c_idx_old in exe_pool.imap(self._split_merge_sampling, args):
            if accepted:
                # if split:
                #    cluster_new = [new cluster 0 induced by init_i, new cluster 1 induced by init_j]
                #    x_split = [[indices of x in new cluster 0], [indices of x in new cluster 1]]
                # if merge:
                #    cluster_new = [new merged cluster]
                #    x_split = [x_restricted]
                cluster_new, x_split = changes

                # Update master <c_launch> and <s_launch>
                for i, cluster in enumerate(cluster_new):   # i = 0, 1 if split; i = 0 if merge
                    c_idx = max(list(c_launch.keys())) + 1
                    c_launch[c_idx] = cluster
                    for idx in x_split[i]:   # if merge, s_changes == [x_restricted]
                        s_launch[idx] = c_idx
                for c_idx in c_idx_old:
                    del c_launch[c_idx]

                # delete previously-rejected information
                for c_idx in c_idx_old:
                    if c_idx in self.rejected_merge_pairs:
                        del self.rejected_merge_pairs[c_idx]

                updated = True
            else:
                if len(c_idx_old) == 2:   # add to the rejected-pair list
                    self.rejected_merge_pairs[c_idx_old[0]].add(c_idx_old[1])
                    self.rejected_merge_pairs[c_idx_old[1]].add(c_idx_old[0])
        exe_pool.close()

        max_post_updated = False

        if updated:   # at least one proposal was accepted
            self.c = copy.deepcopy(c_launch)
            self.s = s_launch.copy()

            logp_post = self.calculate_posterior()
            self.posterior_ncluster.append((logp_post, len(self.c)))
            if self.max_post < logp_post:   # max posterior is updated
                self.max_post = logp_post
                self.max_s = self.s.copy()
                max_post_updated = True
                if self.verbose:
                    print("posterior: %f -> %f" % (self.max_post, logp_post))
        else:
            self.posterior_ncluster.append(self.posterior_ncluster[-1])

        return max_post_updated

    def gibbs_sampling(self):
        # Drop 1 assignment in s
        r = random.randint(0, self.N - 1)
        out_cluster = self.s[r]
        self.c[out_cluster].remove(self.x[r])

        # Delete the cluster from which the data was dropped if it is now empty
        if self.c[out_cluster].ndata == 0:
            del self.c[out_cluster]

        logp_cluster = np.zeros(len(self.c.keys()) + 1)   # cluster assignment probabilities

        # probabilities to be assigned to existing clusters
        for i, cluster in enumerate(self.c.values()):
            if i == out_cluster:
                logp_cluster[i] = np.log(cluster.ndata) - np.log(self.N - 1 + self.alpha) + cluster.logp_x(self.x[r])
            else:
                logp_cluster[i] = np.log(cluster.ndata) - np.log(self.N - 1 + self.alpha) + cluster.logp_x(self.x[r], precomputed=True)

        # probability to be assigned to a new cluster
        logp_cluster[-1] = np.log(self.alpha) - np.log(self.N - 1 + self.alpha) + self.L * self.log_normalizer_const   # TODO: define as constant?

        p_cluster = np.exp(logp_cluster) / np.sum(np.exp(logp_cluster))   # normalized probabilities for assignment
        p_acc_cluster = np.add.accumulate(p_cluster)   # accumulated probabilities

        # assign new cluster based on the p
        rnd = random.random()
        for i, p in enumerate(p_acc_cluster):
            if p >= rnd:
                if i < len(self.c):   # existing cluster
                    c_idx = list(self.c.keys())[i]
                else:   # new cluster
                    c_idx = max(list(self.c.keys())) + 1
                    cluster = Cluster(self.alpha, self.L)
                    self.c[c_idx] = cluster
                break

        self.c[c_idx].add(self.x[r])
        self.s[r] = c_idx

        max_post_updated = False

        if out_cluster != c_idx:   # dropped data is assigned to a different cluster (change in clustering)
            if self.verbose:
                print("%d (%d) -> %d (%d)" % (r, self.s[r], i, c_idx), p_cluster)
                print(self.s)

            # Update parameters and likelihoods of the changed clusters here
            if out_cluster in self.c:   # when not dissapeared
                self.c[out_cluster].calculate_likelihood()
            self.c[c_idx].calculate_likelihood()

            logp_post = self.calculate_posterior()
            self.posterior_ncluster.append((logp_post, len(self.c)))
            if self.max_post < logp_post:
                if self.verbose:
                    print("posterior: %f -> %f #" % (self.max_post, logp_post), self.s)
                self.max_post = logp_post
                self.max_s = self.s.copy()
                max_post_updated = True
        else:
            self.posterior_ncluster.append(self.posterior_ncluster[-1])

        return max_post_updated


def initialize_clustering(vmatrix_fname, alpha=1., n_parallel=1, verbose=False):
    return Clustering(load_vmatrix(vmatrix_fname),
                      alpha,
                      n_parallel,
                      verbose=verbose)


def load_clustering(pkl_fname, n_parallel=1, verbose=False):
    with open(pkl_fname, mode='rb') as f:
        clustering = pickle.load(f)
    clustering.n_parallel = n_parallel
    clustering.verbose = verbose

    print("[INFO] loaded clustering")
    print("[INFO] alpha = %f, N = %d, L = %d, %d iterations so far"
          % (clustering.alpha,
             clustering.N,
             clustering.L,
             len(clustering.posterior_ncluster) - 1))
    print("[INFO] # of clusters = %d (max posterior), %d (current)"
          % (len(set(clustering.max_s)),
             len(clustering.c)))

    return clustering


def run_sampling(clustering,
                 sampling_type,
                 n_iteration,
                 converge_threshold_sm=10000,
                 converge_threshold_gibbs=100000,
                 verbose=False,
                 out_cluster_pickle="clustering.pkl"):

    def output_clustering():
        with open(out_cluster_pickle, 'wb') as f:
            pickle.dump(clustering, f)

    if sampling_type == "both":
        sm_converge = False
        gibbs_converge = False
        no_update = 0

    t_start = time.time()
    t_block_start = time.time()

    for iteration in range(n_iteration):
        if verbose:
            print("## --- step %d --- ##" % iteration)
            print("s =", clustering.s)
            print("c =", list(clustering.c.keys()))

        if ((sampling_type == "sm" or
             (sampling_type == "both" and not sm_converge))):

            # Output temporary pickle files for every 10,000 iterations
            if iteration != 0 and iteration % 10000 == 0:
                t_block_end = time.time()
                output_clustering()
                print("[INFO] Write clustering pickle (SM; %f second for this 10,000 iterations)" % (t_block_end - t_block_start))
                t_block_start = time.time()

            updated = clustering.split_merge_sampling()
            #updated = clustering.split_merge_sampling_parallel()   # multiple proposals (obsolete; too slow)

            if sampling_type == "both":
                # Check convergence
                if not updated:
                    no_update += 1
                    if no_update >= converge_threshold_sm:   # TODO: avoid "pseudo-converge" in the very beggining iterations
                        sm_converge = True
                        no_update = 0
                        print("[INFO] SM converged at %d iterations" % (iteration + 1))   # TODO: output to different pickle files for SM and Gibbs
                        t_block_start = time.time()
                else:
                    no_update = 0

        elif ((sampling_type == "gibbs" or
               (sampling_type == "both" and not gibbs_converge))):

            # Output temporary pickle files for every 100,000 iterations
            if iteration != 0 and iteration % 100000 == 0:
                t_block_end = time.time()
                output_clustering()
                print("[INFO] Write clustering pickle (Gibbs; %f second for this 100,000 iterations)" % (t_block_end - t_block_start))
                t_block_start = time.time()

            updated = clustering.gibbs_sampling()

            if sampling_type == "both":
                # Check convergence
                if not updated:
                    no_update += 1
                    if no_update >= converge_threshold_gibbs:
                        gibbs_converge = True
                        print("[INFO] Gibbs converged at %d iterations" % (iteration + 1))   # TODO: output to different pickle files for SM and Gibbs
                else:
                    no_update = 0

        else:   # both sm and gibbs are converged in the both mode
            break

    # Output final clustering
    output_clustering()

    t_end = time.time()
    print("[INFO] Finished in %f second" % (t_end - t_start))


def main():
    args = load_args()

    # Prepare clustering object
    if args.vmatrix_fname is not None:
        clustering = initialize_clustering(args.vmatrix_fname,
                                           args.alpha,
                                           args.n_parallel,
                                           args.verbose)
    else:
        clustering = load_clustering(args.in_cluster_pickle,
                                     args.n_parallel,
                                     args.verbose)

    # Update the clustering in a specified way
    run_sampling(clustering,
                 args.sampling_type,
                 args.n_iteration,
                 converge_threshold_sm=args.converge_sm,
                 converge_threshold_gibbs=args.converge_gibbs,
                 verbose=args.verbose,
                 out_cluster_pickle=args.out_cluster_pickle)


def load_args():
    parser = argparse.ArgumentParser(
        description=("Dirichlet process mixture model for unit sequence "
                     "clustering."))

    parser.add_argument(
        "--vmatrix_fname",
        default=None,
        help=("Input variant matrix generated by Consed"))

    parser.add_argument(
        "--in_cluster_pickle",
        default="clustering.pkl",
        help=("Input pickle file of a clustering [clustering.pkl]"))

    parser.add_argument(
        "--out_cluster_pickle",
        default="clustering.pkl",
        help=("Output pickle file of a clustering [clustering.pkl]"))

    parser.add_argument(
        "--alpha",
        type=float,
        default=1.,
        help=("Concentration hyperparameter in DPMM [1.]"))

    parser.add_argument(
        "--sampling_type",
        default="sm",
        choices=["sm", "gibbs", "both"],
        help=("Sampling method in DPMM. Split-merge MCMC (sm), Gibbs "
              "sampling (gibbs), or both of them (both). [sm]"))

    parser.add_argument(
        "--n_iteration",
        type=int,
        default=10000,
        help=("Number of iterations in DPMM [10000]"))

    parser.add_argument(
        "--converge_sm",
        type=int,
        default=10000,
        help=("Maximum number of consecuitive iterations with no update on "
              "the posterior probability in Split-merge MCMC. This is the "
              "threshold for convergence. This value depends on inputs. "
              "This is used only in the \"both\" mode. [10000]"))

    parser.add_argument(
        "--converge_gibbs",
        type=int,
        default=100000,
        help=("Maximum number of consecuitive iterations with no update on "
              "the posterior probability in Gibbs sampling. This is the "
              "threshold for convergence. This value depends on inputs. "
              "This is used only in the \"both\" mode. [100000]"))

    parser.add_argument(
        "--n_parallel",
        type=int,
        default=1,
        help=("Number of proposals in each Split-merge MCMC sampling "
              "iteration. The same number of processes will be used. [1]"))

    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help=("Verbose mode [False]"))

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
