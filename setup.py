from setuptools import setup, find_packages

setup(name='fetch_tool',
      version='0.1',
      description='Tool which allows you to fetch RAW read files from the '
                  'European Nucleotide Archive (ENA).',
      url='https://github.com/EBI-Metagenomics/fetch_tool',
      author='Maxim Scheremetjew',
      license='Apache Software License',
      packages=find_packages(exclude=['ez_setup']),
      zip_safe=False,
      python_requires=">=3.4",)
