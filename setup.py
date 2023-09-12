from setuptools import setup, find_packages

""" setup
SUMMARY: setup the package with some infos

EXTENDED SUMMARY: setup define commond variable usefull for the package, as its name, version, packages, and required packages.
"""

setup(
    name="Chemical Pathways Analysis Program",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # List your package's dependencies here
    ],
)