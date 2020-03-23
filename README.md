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


