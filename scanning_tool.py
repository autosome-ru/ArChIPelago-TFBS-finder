#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import sys
import argparse
import subprocess
from sklearn.externals import joblib
from Bio import SeqIO
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


# In[3]:


#### грузим параметры
print("Loading parameters ...")
parser = argparse.ArgumentParser(description='Scanner')
parser.add_argument('-root', action='store', help='Root directory', required=True)
parser.add_argument('-mfa', action='store', help='Multi fasta file', required=True)
parser.add_argument('-f_name', action='store', help='Factor name', required=True)
parser.add_argument('-sarus_h', action='store', help='Sarus PWM scanner home', required=True)
parser.add_argument('-threshold', action='store', help='Result threshold [0-1]', required=True)
parser.add_argument('-step', action='store', help='Step value [1-100]', required=True)
parser.add_argument('-frame', action='store', help='Frame size [200]', required=True)
args = parser.parse_args()


root = args.root #"/home/pavel-kravchenko/basicdir/program_files"
os.chdir(root)
mfa_filename = args.mfa #"sequence_Y.fasta" 
TF_name = args.f_name #"CEBPB" 
threshold = args.threshold
frame = args.frame
step = args.step
pwmdir = root + "/pwmdir_tmp"
featuredir = root + "/featuredir_tmp"
sarus_home = args.sarus_h #"/home/pavel-kravchenko/sarus/releases/sarus-2.0.1.jar"

line = "rm pwmdir_tmp -R; rm featuredir_tmp -R"
p = subprocess.Popen(line, shell=True)
p.wait()

print("Creating directories ...")
if not os.path.exists(pwmdir): # если нет result, то создаем
    os.makedirs(pwmdir)
if not os.path.exists(featuredir): # если нет outputdir, то создаем
    os.makedirs(featuredir)


# In[4]:

print("Processing mfa file ...")
with open("tmp.fasta", "w") as f:
    with open(featuredir + "/headerStr_names.tab", "w") as f_names:
        iterator_fasta = SeqIO.parse(mfa_filename, "fasta")

        for seq in iterator_fasta:
            ind = 0

            if len(seq.seq) > frame:
                for i in range(0,len(seq.seq)-frame,step):
                    tmp_seq = seq.seq[i:i+frame]
                    if all(l not in tmp_seq for l in ["B", "D", "H", "K", "M", "S", "V", "W", "N", "Y"]):
                        f.write(">" + str(seq.id) + "_" + str(ind) + "\n" + str(tmp_seq) + "\n")
                        f_names.write(str(seq.id) + "_" + str(ind) + "\t" + str(tmp_seq) + "\n")
                        ind += 1
            else:
                if all(l not in tmp_seq for l in ["B", "D", "H", "K", "M", "S", "V", "W", "N", "Y"]):
                    f.write(">" + str(seq.id) + "_" + str(ind) + "\n" + str(seq.seq) + "\n")
                    f_names.write(str(seq.id) + "_" + str(ind) + "\t" + str(seq.seq) + "\n")
                    ind += 1


# In[5]:


### извлечь матрицы для фактора в tmp директорию
print("Sarus is scanning ...")
pwm_list = [f for f in [f for f in os.listdir(root) if f.split(".")[-1] == "gz"] if f.split("pwm_files_")[1].split("_HUMAN.tar.gz")[0] == TF_name]
if len(pwm_list) == 0:
    print("This factor is unknown!")
    sys.exit(0)
    
line = "tar xvzf {archive_f} -C {pwmdir}".format(archive_f=root + "/" + pwm_list[0], pwmdir=pwmdir)
p = subprocess.Popen(line, shell=True)
p.wait()

pwm_list = [f for f in os.listdir(pwmdir) if f.split(".")[1] == "pwm"]    
for matrix in pwm_list:
    print(matrix)
    line = 'java -Xmx2G -cp {sarus_home} ru.autosome.SARUS {mfa_file} {matrix} --skipn --show-non-matching --output-scoring-mode score besthit | grep -v \> > {out_file}'.format(sarus_home = sarus_home,
    mfa_file = "tmp.fasta",
    matrix = '/'.join([pwmdir, matrix]),
    out_file = '/'.join([featuredir, matrix.split(".")[0] + ".tab"]))
    print(line)
    p = subprocess.Popen(line, shell=True)
    p.wait()


# In[6]:


feature_list = [f for f in os.listdir(featuredir) if f.split(".")[1] == "tab" and f != "headerStr_names.tab"]   
string = "paste <( cut -f 1-2 {out_tab_file} )".format(out_tab_file=featuredir + "/headerStr_names.tab")
for fn in feature_list:
    #print(fn)
    string += ' <( cut -f 1 {file_name} )'.format(file_name=featuredir + "/" + fn)
string += ' > {featuredir}/out_sure.csv'.format(featuredir=featuredir)
#print(string)

with open("code.sh", "w") as f:
    f.write(string)

line = 'bash code.sh'
p = subprocess.Popen(line, shell=True)
p.wait()

df = pd.read_csv(featuredir + "/" + "out_sure.csv", header=None, sep='\t')
#print(df.head())

df1 = df.iloc[:,2:]
df1.columns = [x for x in range(4,24)]
#print(df1)
 
model_filename = [f for f in [f for f in os.listdir(root) if os.path.splitext(f)[1] == '.sav'] if f.split("finalized_model_")[1].split("_")[0] == TF_name] 
#print(archive_f)

print("Predicting ...")
MODEL = joblib.load(root + "/" + model_filename[0])
MODEL

Y_test_predicted_proba = MODEL.predict_proba(df1)[:, 1]
#print(Y_test_predicted_proba[:10])

df["Y_test_predicted_proba"] = Y_test_predicted_proba
#out_df = df[df["Y_test_predicted_proba"] > 0.6].sort_values(by=["Y_test_predicted_proba"], ascending=False)
out_df = df[df["Y_test_predicted_proba"] > threshold]


# In[8]:

print("Saving results into out_df.csv file ...")
pd.DataFrame.to_csv(out_df, "out_df.csv",  sep='\t', header=True)


print("Done!")

