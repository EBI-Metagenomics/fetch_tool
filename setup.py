import os
import sys

from setuptools import find_packages, setup

_base = os.path.dirname(os.path.abspath(__file__))
_requirements = os.path.join(_base, "requirements.txt")
_requirements_test = os.path.join(_base, "requirements-test.txt")

version = "0.8.0"

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
    description="Utility to fetch public and private RAW read and assembly "
    "files from the ENA",
    author="Microbiome Informatics Team",
    url="https://github.com/EBI-Metagenomics/fetch_tool",
    packages=find_packages(),
    install_requires=install_requirements,
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "fetch-assembly-tool=fetchtool.fetch_assemblies:main",
            "fetch-read-tool=fetchtool.fetch_reads:main",
        ]
    },
    tests_require=test_requirements,
    test_suite="tests",
    setup_requires=["pytest-runner"],
)
