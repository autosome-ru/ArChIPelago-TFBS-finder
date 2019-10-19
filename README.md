# TFBS-finder

This is a repo of transcription factors binding sites prediction project

This script provides a TFBS search based on a PWMs' unified classifier prediction. 

## Before you start

Make sure that you have installed all components:
<ul>
<li>Python 3.6 or upper https://www.python.org/ + (subprocess, matplotlib,pandas,numpy,itertools,seaborn,sklearn,pickle,csv,xgboost,lightgbm ,catboost,collections,operator)
<li>Biopython 1.70 or upper http://biopython.org/
<li>SPRY-SARUS http://autosome.ru/ChIPMunk/
</ul>


## Getting started

### Installation

First of all you have to ```clone``` this directory</br></br>
```git clone https://github.com/Pavel-Kravchenko/TFBS-finder/```</br></br>
Then ```cd``` in TFBS-finder</br></br>
```cd TFBS-finder```</br></br>

Now you are ready to start.
Run the script with your .mfa file for selected TF. 
``` python scanning_tool.py -mfa [your_mfa_file] -f_name [factor_of_interest] -threshold 0.5 -step 50```


## Contact me

Feel free to contact me for any suggestions or critique.

Email: pavel-kravchenk0@yandex.ru 

Site: http://kodomo.fbb.msu.ru/~pavel-kravchenko/index.html 
