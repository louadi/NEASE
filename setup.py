from setuptools import setup, find_packages


# The text of the README file
README = (HERE / "README.md").read_text()


setup(name='nease',
      version='0.0.1',
      description='NEASE: A network-based approach for the enrichment of alternative splicing events',
      long_description=README,
      url='https://github.com/louadi/NEASE',
      author='Zakaria Louadi',
      author_email='Zakaria.louadi@tum.de',
      license='GPLv3',
      packages=find_packages(),
      python_requires='>=3.6',
      install_requires=[
        'pandas>=1.2.1',
        'numpy>=1.20.2',
        'networkx=2.4',
        'plotly>=4.14.3',
        'scipy>=1.6.2',
          
        'statsmodels>=0.9.15',
        'seaborn>=0.9.0',
        'mygene>=3.1.0',
        'scikit_learn>=0.22',
          
    statsmodels==0.12.2
        
    ],

)