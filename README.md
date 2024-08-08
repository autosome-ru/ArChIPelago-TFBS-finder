# TFBS-finder: ArChIPelago — Aggregatin of multiple position weight matrices and ChIP-seq with machinE LeArninG for predictiOn of transcription factors binding sites


This is a repo of transcription factors binding sites prediction project.

This script provides a TFBS search with the help of unified classifier built on mononucleotide PWMs. 

## Before you start

Make sure that you have installed all components:
<ul>
<li>Python 3.6 or upper https://www.python.org/ + (subprocess, matplotlib,pandas,numpy,itertools,seaborn,sklearn,pickle,csv,xgboost,lightgbm ,catboost,collections,operator)
<li>Biopython 1.70 or upper http://biopython.org/
<li>SPRY-SARUS http://autosome.ru/ChIPMunk/
</ul>


## Getting started

### Installation

First of all you have to ```clone``` this directory.
</br></br>
```git clone https://github.com/Pavel-Kravchenko/TFBS-finder/```
</br></br>
Then ```cd``` in TFBS-finder.
</br></br>
```cd TFBS-finder```
</br></br>

Now you are ready to start.
Run the script with your .mfa file for selected TF. 
</br></br>
``` python scanning_tool.py -mfa [your_mfa_file] -f_name [factor_of_interest] -threshold [0-1]```
</br></br>
</br></br>
Example1: ```python scanning_tool.py -mfa random_10000.fasta -tf_name all -threshold 0.5 -pwm m -model RandomForestClassifier```
</br></br>
Example2: ```python scanning_tool.py -mfa random_10000.fasta -tf_name CTCF -threshold 0.5 -pwm d -model LogisticRegression```
</br></br>
Example3: ```python scanning_tool.py -mfa random_10000.fasta -tf_name ANDR -threshold 0.5 -pwm md -model BaggingClassifier_LogisticRegression```
</br></br>

## Contact me

Feel free to contact me for any suggestions or additional information on the project.

Email: pavel-kravchenk0[@]yandex[dot]ru 

Site: http://kodomo.fbb.msu.ru/~pavel-kravchenko/index.html 
