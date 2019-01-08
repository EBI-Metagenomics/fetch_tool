from setuptools import setup, find_packages
import os
import sys

_base = os.path.dirname(os.path.abspath(__file__))
_requirements = os.path.join(_base, 'requirements.txt')
_requirements_test = os.path.join(_base, 'requirements-test.txt')

version = '0.1.0'

install_requirements = []
with open(_requirements) as f:
    install_requirements = f.read().splitlines()

test_requirements = []
if "test" in sys.argv:
    with open(_requirements_test) as f:
        test_requirements = f.read().splitlines()

setup(name='fetch-tool',
      version=version,
      description='Utility to fetch public and private RAW read and assembly '
                  'files from the ENA',
      author='Maxim Scheremetjew',
      url='https://github.com/EBI-Metagenomics/fetch_tool',
      packages=find_packages(),
      install_requires=install_requirements,
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'fetch-assembly-tool=src.fetch_assemblies:main',
              'fetch-read-tool=src.fetch_reads:main'
          ]
      },
      tests_require=test_requirements,
      test_suite="tests",
      setup_requires=['pytest-runner'],
      )
