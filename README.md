Mummichog 3 dev
===============

Mummichog is a Python program for analyzing data from high throughput, untargeted metabolomics.
It leverages the organization of metabolic networks to predict functional activity directly from feature tables,
bypassing metabolite identification. The version 2 is hosted at:

https://github.com/shuzhao-li/mummichog

This is version 3 under development.


*To run test:*

shuzhao@canyon:~/li.github/mummichog3$ python3 -m mummichog.main -f mummichog/tests/testdata0710.txt -o t3


## Separating out of core package 

* Metabolic models (via azimuth/)

* Visualization (via report/)

* Annotation is optional


## The mummichog suite include

* mummichog(3): core algorithm package for pathway/network analysis

* cloud-mummichog: server and worker (RESTful) implementations

* Azimuth DB: the chemical database for biology, including metabolic models

* metDataModel: data models for metabolomics, used by mummichog and Azimuth DB

* mass2chem: common utilities in interpreting mass spectrometry data, annotation

* massBrowser: visualization using js


## Dev notes

metDataModel is renamed from azimuth-metabolomics.

mummichog includes a default metabolic model, but tries to connect to Azimuth DB for latest models.

Message broker is in cloud-mummichog.
