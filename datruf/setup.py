import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="datruf",
    version="0.1.0",
    description=("DAzzler Tandem Repeat Unit Finding module"),
    long_description=read("README.md"),
    author="Yoshihiko Suzuki",
    author_email="ys.neoteny@gmail.com",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords='datruf',
    install_requires=[
        "numpy",
        "pandas"
    ],
    packages=[
        "datruf"
    ]
)
