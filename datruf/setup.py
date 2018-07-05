from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('requirements.txt') as f:
    requirements = f.read().strip().split('\n')


setup(
    name="datruf",
    version="0.1.0",
    description=("DAzzler Tandem Repeat Unit Finding module"),
    long_description=readme,
    author="Yoshihiko Suzuki",
    author_email="ys.neoteny@gmail.com",
    install_requires=requirements,
    packages=find_packages(exclude=('tests', 'docs'))
)
