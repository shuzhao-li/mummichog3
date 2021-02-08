Mummichog 3 dev
===============

Mummichog is a Python program for analyzing data from high throughput, untargeted metabolomics.
It leverages the organization of metabolic networks to predict functional activity directly from feature tables,
bypassing metabolite identification. The version 2 is hosted at:

https://github.com/shuzhao-li/mummichog

This is version 3 under development.

## set up env for development (Python3, using virtualenv on Linux)

sudo apt install python3-dev python3-venv

sudo pip3 install virtualenv

virtualenv env

source env/bin/activate

A few libraries used for mummichog, for example:

(env) $ pip install scipy matplotlib xlsxwriter

(env) $ pip install networkx==1.10

(env) $ deactivate

*To run test:*

shuzhao@canyon:~/li.github/mummichog3$ python3 -m mummichog.main -f mummichog/tests/testdata0710.txt -o t3


## Separating out of core package 

* Metabolic models (via azimuth/)

* Visualization (via report/)

* Annotation is optional


## The mummichog suite includes

* mummichog(3): core algorithm package for pathway/network analysis

* cloud-mummichog: server and worker (RESTful) implementations

* Azimuth DB: the chemical database for biology, including metabolic models

* metDataModel: data models for metabolomics, used by mummichog and Azimuth DB

* mass2chem: common utilities in interpreting mass spectrometry data, annotation

* massBrowser: visualization using js


## Dev notes

Be aware that mummichog (> 3.0.3) imports models from metDataModel and adducts from mass2chem,

mummichog includes a default metabolic model, but tries to connect to Azimuth DB for latest models.

Message broker is in cloud-mummichog. Francisco is working on similar functions in mummichog3-api (JAX BitBucket).
