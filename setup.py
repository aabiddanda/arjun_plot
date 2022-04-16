"""Setup for the arjun_plot library."""

from setuptools import setup

version = "0.0.1"

required = open("requirements.txt").read().split("\n")

setup(
    name="arjun_plot",
    version=version,
    description="Custom matplotlib plotting routines",
    author="aabiddanda",
    author_email="aabiddanda@gmail.com",
    url="https://github.com/aabiddanda/arjun_plot",
    packages=["arjun_plot"],
    install_requires=required,
    long_description="See " + "https://github.com/aabiddanda/arjun_plot",
    license="MIT",
)
