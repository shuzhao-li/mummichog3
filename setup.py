from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
  name='mummichog',
  version='3.0.4',

  author='Shuzhao Li, Francisco Castellanos, Andrei Todor',
  author_email='shuzhao.li@gmail.com',
  description='Pathway and network analysis for metabolomics data',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='https://github.com/shuzhao-li/mummichog3',
  license='BSD 3-Clause',
  keywords='metabolomics analysis bioinformatics mass spectrometry systems biology',

  classifiers=[
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=find_packages(),
  data_files=[ ('mummichog/tests', ['mummichog/tests/testdata0710.txt']) ],
  include_package_data=True,
  zip_safe=True,
  entry_points = {
        'console_scripts': ['mummichog=mummichog.command_line:main'],
    },

  python_requires='>=3.4',
  install_requires=[
    'metDataModel',
    'mass2chem',
    'matplotlib',
    'networkx>=1,<2',
    'numpy',
    'scipy',
    'xlsxwriter',
  ],

)
