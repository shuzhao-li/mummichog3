from setuptools import setup, find_packages

with open("mummichog/__init__.py") as f:
    exec([x for x in f.readlines() if '__version__' in x][0])

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
  name='mummichog',
  version=__version__,

  author='Shuzhao Li',
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
    'networkx',
    'numpy',
    'scipy',
    'xlsxwriter',
  ],

)
