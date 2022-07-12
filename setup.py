import os
import sys

from setuptools import find_packages, setup

_base = os.path.dirname(os.path.abspath(__file__))
_requirements = os.path.join(_base, "requirements.txt")
_requirements_test = os.path.join(_base, "requirements-test.txt")

version = "0.8.1"

install_requirements = []
with open(_requirements) as f:
    install_requirements = f.read().splitlines()

test_requirements = []
if "test" in sys.argv:
    with open(_requirements_test) as f:
        test_requirements = f.read().splitlines()

setup(
    name="fetch-tool",
    version=version,
    description="Microbiome Informatics ENA Fetch Tool",
    long_description="Utility to fetch public and private RAW read and assembly files from the ENA",
    author="Microbiome Informatics Team",
    url="https://github.com/EBI-Metagenomics/fetch_tool",
    packages=find_packages(),
    install_requires=install_requirements,
    include_package_data=True,
    python_requires=">=3.7",
    license="Apache Software License",
    tests_require=test_requirements,
    test_suite="tests",
    setup_requires=["pytest-runner"],
    entry_points={
        "console_scripts": [
            "fetch-assembly-tool=fetchtool.fetch_assemblies:main",
            "fetch-read-tool=fetchtool.fetch_reads:main",
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
)
