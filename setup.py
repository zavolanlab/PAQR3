"""PAQR3 package definition"""

from setuptools import (setup, find_packages)

from paqr3.version import __version__

# # Read long description from file
# with open("README.md", "r", encoding="utf-8") as fh:
#     LONG_DESCRIPTION = fh.read()

setup(
    name="paqr3",
    version=__version__,
    description=(
        "Quantification of poly(A) sites on standard RNA-Seq data."
    ),
    url="https://github.com/zavolanlab/paqr3",
    author="Máté Balajti",
    author_email="mate.balajti@unibas.ch",
    maintainer="Máté Balajti",
    maintainer_email="mate.balajti@unibas.ch",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
    ],
    entry_points={
        'console_scripts': [
            'paqr3 = paqr3.cli:main',
        ],
    },
    keywords=[
        'bioinformatics',
        'polya',
        'quantification',
    ],
    project_urls={
        "Repository": "https://github.com/zavolanlab/paqr3",
        "Tracker": "https://github.com/zavolanlab/paqr3/issues",
    },
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[
        "setuptools_git == 1.2",
    ],
)