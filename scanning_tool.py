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
<<<<<<< HEAD
matplotlib.use('agg')
=======
>>>>>>> 51296126fcc189628509f4fc3e053fc75b3d4dbf
import matplotlib.pyplot as plt


# In[3]:


#### грузим параметры
print("Loading parameters ...")
parser = argparse.ArgumentParser(description='Scanner')
parser.add_argument('-root', action='store', help='Root directory', required=False, default=".")
parser.add_argument('-mfa', action='store', help='Multi fasta file', required=True)
<<<<<<< HEAD
parser.add_argument('-tf_name', action='store', help='Factor name', required=True)
parser.add_argument('-sarus_h', action='store', help='Sarus PWM scanner home', required=False, default="./sarus/releases/sarus-2.0.1.jar")
parser.add_argument('-threshold', action='store', help='Result threshold [0-1]', required=True)
parser.add_argument('-step', action='store', help='Step value (def. 10)', required=False, default=10)
parser.add_argument('-frame', action='store', help='Frame size (def. 200)', required=False, default=200)
parser.add_argument('-cov_plt', action='store', help='Do you want to see a coverage plot [1/0] (def. 0)', required=False, default=0)
parser.add_argument('-pwm', action='store', help='Do you want to use mononucleotide, dinucleotide or both pwm [m/d/md] (def. m)', required=False, default="m")
parser.add_argument('-tf_av', action='store', help='Print all available TF [y/n] (def. n)', required=False, default="n")
=======
parser.add_argument('-f_name', action='store', help='Factor name', required=True)
parser.add_argument('-sarus_h', action='store', help='Sarus PWM scanner home', required=False, default="./sarus/releases/sarus-2.0.1.jar")
parser.add_argument('-threshold', action='store', help='Result threshold [0-1]', required=True)
parser.add_argument('-step', action='store', help='Step value [1-100]', required=False, default=100)
parser.add_argument('-frame', action='store', help='Frame size [200]', required=False, default=200)
>>>>>>> 51296126fcc189628509f4fc3e053fc75b3d4dbf
args = parser.parse_args()


root = args.root #"/home/pavel-kravchenko/basicdir/program_files"
os.chdir(root)
<<<<<<< HEAD

mfa_filename = args.mfa #"sequence_Y.fasta" 
TF_name = args.tf_name #"CEBPB" 
coverage_plot = int(args.cov_plt)
threshold = float(args.threshold)
frame = int(args.frame)
step = int(args.step)
pwm = args.pwm
basicdir = root + "/Factors/" + TF_name
pwmdir = root + "/Factors/" + TF_name + "/pwmdir_tmp"
dpwmdir = root + "/Factors/" + TF_name + "/dpwmdir_tmp"
featuredir = root + "/Factors/" + TF_name + "/featuredir_tmp"
sarus_home = args.sarus_h #"/home/pavel-kravchenko/sarus/releases/sarus-2.0.1.jar"
results = root + "/Results"

if args.tf_av == "y":
    print("TF available ...")
    for i in ['ANDR_HUMAN',
	 'CEBPB_HUMAN',
	 'ESR1_HUMAN',
	 'GATA1_HUMAN',
	 'MYC_HUMAN',
	 'PO5F1_HUMAN',
	 'REST_HUMAN',
	 'RXRA_HUMAN',
	 'SPI1_HUMAN',
	 'STAT1_HUMAN',
	 'CEBPA_HUMAN',
	 'CTCF_HUMAN',
	 'FOXA1_HUMAN',
	 'GCR_HUMAN',
	 'P53_HUMAN',
	 'PPARG_HUMAN',
	 'RUNX1_HUMAN',
	 'SOX2_HUMAN',
	 'STA5A_HUMAN',
	 'TF65_HUMAN']:
        print(i)
    print(" ")

=======
mfa_filename = args.mfa #"sequence_Y.fasta" 
TF_name = args.f_name #"CEBPB" 
threshold = float(args.threshold)
frame = int(args.frame)
step = int(args.step)
pwmdir = root + "/pwmdir_tmp"
featuredir = root + "/featuredir_tmp"
sarus_home = args.sarus_h #"/home/pavel-kravchenko/sarus/releases/sarus-2.0.1.jar"

line = "rm pwmdir_tmp -R; rm featuredir_tmp -R"
p = subprocess.Popen(line, shell=True)
p.wait()
>>>>>>> 51296126fcc189628509f4fc3e053fc75b3d4dbf

print("Creating directories ...")
if not os.path.exists(pwmdir): # если нет result, то создаем
    os.makedirs(pwmdir)
<<<<<<< HEAD
if not os.path.exists(dpwmdir): # если нет result, то создаем
    os.makedirs(dpwmdir)
if not os.path.exists(featuredir): # если нет outputdir, то создаем
    os.makedirs(featuredir)
if not os.path.exists(results): # если нет outputdir, то создаем
    os.makedirs(results)

=======
if not os.path.exists(featuredir): # если нет outputdir, то создаем
    os.makedirs(featuredir)
>>>>>>> 51296126fcc189628509f4fc3e053fc75b3d4dbf


# In[4]:

print("Processing mfa file ...")
with open("tmp.fasta", "w") as f:
    with open(featuredir + "/headerStr_names.tab", "w") as f_names:
<<<<<<< HEAD
        with open(featuredir + "/headerStr_names_init.tab", "w") as f_names_init:
            iterator_fasta = SeqIO.parse(mfa_filename, "fasta")
            for seq in iterator_fasta:
                ind = 0
                f_names_init.write(str(seq.id) + "\t" + str(seq.seq) + "\n")
                if len(seq.seq) > frame:
                    for i in range(0,len(seq.seq)-frame,step):
                        tmp_seq = seq.seq[i:i+frame]
                        if all(l not in tmp_seq for l in ["B", "D", "H", "K", "M", "S", "V", "W", "N", "Y"]):
                            f.write(">" + str(seq.id) + "@" + str(ind) + "\n" + str(tmp_seq) + "\n")
                            f_names.write(str(seq.id) + "@" + str(ind) + "\t" + str(tmp_seq) + "\n")
                            ind += 1
                else:
                    if all(l not in tmp_seq for l in ["B", "D", "H", "K", "M", "S", "V", "W", "N", "Y"]):
                        f.write(">" + str(seq.id) + "@" + str(ind) + "\n" + str(seq.seq) + "\n")
                        f_names.write(str(seq.id) + "@" + str(ind) + "\t" + str(seq.seq) + "\n")
                        ind += 1

# In[5]:
### извлечь матрицы для фактора в tmp директорию

print("Sarus is scanning ...")
if pwm == "m":
    pwm_list = [f for f in [f for f in os.listdir(basicdir) if f.split(".")[-1] == "gz"] if f[0] == "p"]
    if len(pwm_list) == 0:
        print("This factor is unknown!")
        sys.exit(0)

    if len(os.listdir(pwmdir)) == 0:
	    line = "tar xvzf {archive_f} -C {pwmdir}".format(archive_f=basicdir + "/" + pwm_list[0], pwmdir=pwmdir)
	    p = subprocess.Popen(line, shell=True)
	    p.wait()

    pwm_list = [f for f in os.listdir(pwmdir) if f.split(".")[1] == "pwm"]    
    for matrix in pwm_list:
        print(matrix)
        line = 'java -Xmx2G -cp {sarus_home} ru.autosome.SARUS {mfa_file} {matrix} --skipn --show-non-matching --output-scoring-mode score besthit | grep -v \> > {out_file}'.format(sarus_home = sarus_home,
        mfa_file = "tmp.fasta",
        matrix = '/'.join([pwmdir, matrix]),
        out_file = '/'.join([featuredir, matrix.split(".")[0] + "_mono.tab"]))
        print(line)
        p = subprocess.Popen(line, shell=True)
        p.wait()


    # In[6]:
    feature_list = [f for f in os.listdir(featuredir) if f.split("_mono.")[-1] == "tab" and f != "headerStr_names.tab" and f != "headerStr_names_init.tab"]   
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


    df = pd.read_csv(featuredir + "/out_sure.csv", header=None, sep='\t')
    #print(df.head())

    df1 = df.iloc[:,2:]
    print(df1)
    df1.columns = [x for x in range(4,24)]
    #print(df1)

    model_filename = [f for f in [f for f in os.listdir(basicdir) if f.split("_mono.")[-1] == 'sav'] if f.split("XGBClassifier")[-1] == "_mono.sav"] 
    #print(archive_f)

    print("Predicting ...")
    MODEL = joblib.load(basicdir + "/" + model_filename[0])
    MODEL

    Y_test_predicted_proba = MODEL.predict_proba(df1)[:, 1]
    #print(Y_test_predicted_proba[:10])

    df["Y_test_predicted_proba"] = Y_test_predicted_proba
    #out_df = df[df["Y_test_predicted_proba"] > 0.6].sort_values(by=["Y_test_predicted_proba"], ascending=False)
    out_df = df[df["Y_test_predicted_proba"] > threshold]


    if coverage_plot == 1:
        data = {}
        for i in out_df[0]:
            name = i.split("@")[0]
            ind = i.split("@")[-1]
            #print(i)
            if name in data.keys():
                data[name].append([float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind])
            else:
                data[name] = [[float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind]]
            #print(data)

        for i, name in enumerate(list(data.keys())):
            #print(i, name)
            y, x = zip(*data[name])
            #print(x)
            #print(y)
            plt.plot(x, y, 'o')
            plt.xticks(rotation=90)
            plt.title(name)
            plt.xlabel('position')
            plt.ylabel('score')
            plt.tight_layout()
            plt.savefig(results + "/coverage_plot_TF_{TF_name}_threshold_{threshold}_seq_{i}_mono.pdf".format(TF_name=TF_name, threshold=threshold, i=i), dpi=100)
            plt.close()


    # In[8]:

    print("Saving results into out_df_to_cols.csv file ...")
    pd.DataFrame.to_csv(out_df, results + "/out_df_to_cols_TF_{TF_name}_threshold_{threshold}_mono.csv".format(TF_name=TF_name, threshold=threshold),  sep='\t', header=True)

    print("Saving results into out_df_to_rows.csv file ...")

    data = {}
    for i in out_df[0]:
        name = i.split("@")[0]
        ind = i.split("@")[-1]
        #print(i)
        if name in data.keys():
            data[name].append([float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind])
        else:
            data[name] = [[float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind]]
    #print(data)

    rows = []
    cols = []
    for name in list(data.keys()):
        #print(i, name)
        y, x = zip(*data[name])
        rows.append(y)
        for i in x:
            if i not in cols:
                cols.append(i)

    df3 = pd.DataFrame(rows)
    df3.columns = cols
    #df3

    df_init = pd.read_csv(featuredir + "/headerStr_names_init.tab", header=None, sep="\t")
    df_init.columns = ["Name", "Seq"]
    #df_init

    res = pd.concat([df_init, df3], axis=1, sort=False)
    res
    pd.DataFrame.to_csv(res, results + "/out_df_to_rows_TF_{TF_name}_threshold_{threshold}_mono.csv".format(TF_name=TF_name, threshold=threshold),  sep='\t', header=True)
    print("Done!")

if pwm == "d":
    pwm_list = [f for f in [f for f in os.listdir(basicdir) if f.split(".")[-1] == "gz"] if f[0] == "d"]
    print(pwm_list)
    if len(pwm_list) == 0:
        print("This factor is unknown!")
        sys.exit(0)

    if len(os.listdir(dpwmdir)) == 0:
	    line = "tar xvzf {archive_f} -C {dpwmdir}".format(archive_f=basicdir + "/" + pwm_list[0], dpwmdir=dpwmdir)
	    p = subprocess.Popen(line, shell=True)
	    p.wait()

    pwm_list = [f for f in os.listdir(dpwmdir) if f.split(".")[1] == "dpwm"]    
    for matrix in pwm_list:
        print(matrix)
        line = 'java -Xmx2G -cp {sarus_home} ru.autosome.di.SARUS {mfa_file} {matrix} --skipn --show-non-matching --output-scoring-mode score besthit | grep -v \> > {out_file}'.format(sarus_home = sarus_home,
        mfa_file = "tmp.fasta",
        matrix = '/'.join([dpwmdir, matrix]),
        out_file = '/'.join([featuredir, matrix.split(".")[0] + "_di.tab"]))
        print(line)
        p = subprocess.Popen(line, shell=True)
        p.wait()


    # In[6]:
    feature_list = [f for f in os.listdir(featuredir) if f.split("_di.")[-1] == "tab" and f != "headerStr_names.tab" and f != "headerStr_names_init.tab"]   
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


    df = pd.read_csv(featuredir + "/out_sure.csv", header=None, sep='\t')
    #print(df.head())

    df1 = df.iloc[:,2:]
    print(df1)
    df1.columns = [x for x in range(4,24)]
    #print(df1)

    model_filename = [f for f in [f for f in os.listdir(basicdir) if f.split("_di.")[-1] == 'sav'] if f.split("XGBClassifier")[-1] == "_di.sav"] 
    #print(archive_f)

    print("Predicting ...")
    MODEL = joblib.load(basicdir + "/" + model_filename[0])
    MODEL

    Y_test_predicted_proba = MODEL.predict_proba(df1)[:, 1]
    #print(Y_test_predicted_proba[:10])

    df["Y_test_predicted_proba"] = Y_test_predicted_proba
    #out_df = df[df["Y_test_predicted_proba"] > 0.6].sort_values(by=["Y_test_predicted_proba"], ascending=False)
    out_df = df[df["Y_test_predicted_proba"] > threshold]


    if coverage_plot == 1:
        data = {}
        for i in out_df[0]:
            name = i.split("@")[0]
            ind = i.split("@")[-1]
            #print(i)
            if name in data.keys():
                data[name].append([float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind])
            else:
                data[name] = [[float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind]]
            #print(data)

        for i, name in enumerate(list(data.keys())):
            #print(i, name)
            y, x = zip(*data[name])
            #print(x)
            #print(y)
            plt.plot(x, y, 'o')
            plt.xticks(rotation=90)
            plt.title(name)
            plt.xlabel('position')
            plt.ylabel('score')
            plt.tight_layout()
            plt.savefig(results + "/coverage_plot_TF_{TF_name}_threshold_{threshold}_seq_{i}_di.pdf".format(TF_name=TF_name, threshold=threshold, i=i), dpi=100)
            plt.close()


    # In[8]:

    print("Saving results into out_df_to_cols.csv file ...")
    pd.DataFrame.to_csv(out_df, results + "/out_df_to_cols_TF_{TF_name}_threshold_{threshold}_di.csv".format(TF_name=TF_name, threshold=threshold),  sep='\t', header=True)

    print("Saving results into out_df_to_rows.csv file ...")

    data = {}
    for i in out_df[0]:
        name = i.split("@")[0]
        ind = i.split("@")[-1]
        #print(i)
        if name in data.keys():
            data[name].append([float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind])
        else:
            data[name] = [[float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind]]
    #print(data)

    rows = []
    cols = []
    for name in list(data.keys()):
        #print(i, name)
        y, x = zip(*data[name])
        rows.append(y)
        for i in x:
            if i not in cols:
                cols.append(i)

    df3 = pd.DataFrame(rows)
    df3.columns = cols
    #df3

    df_init = pd.read_csv(featuredir + "/headerStr_names_init.tab", header=None, sep="\t")
    df_init.columns = ["Name", "Seq"]
    #df_init

    res = pd.concat([df_init, df3], axis=1, sort=False)
    res
    pd.DataFrame.to_csv(res, results + "/out_df_to_rows_TF_{TF_name}_threshold_{threshold}_di.csv".format(TF_name=TF_name, threshold=threshold),  sep='\t', header=True)
    print("Done!")



if pwm == "md":
    pwm_list = [f for f in [f for f in os.listdir(basicdir) if f.split(".")[-1] == "gz"] if f[0] == "d"]
    print(pwm_list)
    if len(pwm_list) == 0:
        print("This factor is unknown!")
        sys.exit(0)

    if len(os.listdir(dpwmdir)) == 0:
	    line = "tar xvzf {archive_f} -C {dpwmdir}".format(archive_f=basicdir + "/" + pwm_list[0], dpwmdir=dpwmdir)
	    p = subprocess.Popen(line, shell=True)
	    p.wait()

    pwm_list = [f for f in os.listdir(dpwmdir) if f.split(".")[1] == "dpwm"]    
    for matrix in pwm_list:
        print(matrix)
        line = 'java -Xmx2G -cp {sarus_home} ru.autosome.di.SARUS {mfa_file} {matrix} --skipn --show-non-matching --output-scoring-mode score besthit | grep -v \> > {out_file}'.format(sarus_home = sarus_home,
        mfa_file = "tmp.fasta",
        matrix = '/'.join([dpwmdir, matrix]),
        out_file = '/'.join([featuredir, matrix.split(".")[0] + "_di.tab"]))
        print(line)
        p = subprocess.Popen(line, shell=True)
        p.wait()

    pwm_list = [f for f in [f for f in os.listdir(basicdir) if f.split(".")[-1] == "gz"] if f[0] == "p"]
    if len(pwm_list) == 0:
        print("This factor is unknown!")
        sys.exit(0)

    if len(os.listdir(pwmdir)) == 0:
	    line = "tar xvzf {archive_f} -C {pwmdir}".format(archive_f=basicdir + "/" + pwm_list[0], pwmdir=pwmdir)
	    p = subprocess.Popen(line, shell=True)
	    p.wait()

    pwm_list = [f for f in os.listdir(pwmdir) if f.split(".")[1] == "pwm"]    
    for matrix in pwm_list:
        print(matrix)
        line = 'java -Xmx2G -cp {sarus_home} ru.autosome.SARUS {mfa_file} {matrix} --skipn --show-non-matching --output-scoring-mode score besthit | grep -v \> > {out_file}'.format(sarus_home = sarus_home,
        mfa_file = "tmp.fasta",
        matrix = '/'.join([pwmdir, matrix]),
        out_file = '/'.join([featuredir, matrix.split(".")[0] + "_mono.tab"]))
        print(line)
        p = subprocess.Popen(line, shell=True)
        p.wait()


    # In[6]:
    feature_list = [f for f in os.listdir(featuredir) if f.split("_di.")[-1] == "tab" and f != "headerStr_names.tab" and f != "headerStr_names_init.tab"]   
    string = "paste <( cut -f 1-2 {out_tab_file} )".format(out_tab_file=featuredir + "/headerStr_names.tab")
    for fn in feature_list:
        #print(fn)
        string += ' <( cut -f 1 {file_name} )'.format(file_name=featuredir + "/" + fn)
    
    feature_list = [f for f in os.listdir(featuredir) if f.split("_mono.")[-1] == "tab" and f != "headerStr_names.tab" and f != "headerStr_names_init.tab"]   
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


    df = pd.read_csv(featuredir + "/out_sure.csv", header=None, sep='\t')
    #print(df.head())

    df1 = df.iloc[:,2:]
    print(df1)
    df1.columns = [x for x in range(4,44)]
    #print(df1)

    model_filename = [f for f in [f for f in os.listdir(basicdir) if f.split("_di.")[-1] == 'sav'] if f.split("XGBClassifier")[-1] == "_mono_di.sav"] 
    #print(archive_f)

    print("Predicting ...")
    MODEL = joblib.load(basicdir + "/" + model_filename[0])
    MODEL

    Y_test_predicted_proba = MODEL.predict_proba(df1)[:, 1]
    #print(Y_test_predicted_proba[:10])

    df["Y_test_predicted_proba"] = Y_test_predicted_proba
    #out_df = df[df["Y_test_predicted_proba"] > 0.6].sort_values(by=["Y_test_predicted_proba"], ascending=False)
    out_df = df[df["Y_test_predicted_proba"] > threshold]


    if coverage_plot == 1:
        data = {}
        for i in out_df[0]:
            name = i.split("@")[0]
            ind = i.split("@")[-1]
            #print(i)
            if name in data.keys():
                data[name].append([float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind])
            else:
                data[name] = [[float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind]]
            #print(data)

        for i, name in enumerate(list(data.keys())):
            #print(i, name)
            y, x = zip(*data[name])
            #print(x)
            #print(y)
            plt.plot(x, y, 'o')
            plt.xticks(rotation=90)
            plt.title(name)
            plt.xlabel('position')
            plt.ylabel('score')
            plt.tight_layout()
            plt.savefig(results + "/coverage_plot_TF_{TF_name}_threshold_{threshold}_seq_{i}_mono_di.pdf".format(TF_name=TF_name, threshold=threshold, i=i), dpi=100)
            plt.close()


    # In[8]:

    print("Saving results into out_df_to_cols.csv file ...")
    pd.DataFrame.to_csv(out_df, results + "/out_df_to_cols_TF_{TF_name}_threshold_{threshold}_mono_di.csv".format(TF_name=TF_name, threshold=threshold),  sep='\t', header=True)

    print("Saving results into out_df_to_rows.csv file ...")

    data = {}
    for i in out_df[0]:
        name = i.split("@")[0]
        ind = i.split("@")[-1]
        #print(i)
        if name in data.keys():
            data[name].append([float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind])
        else:
            data[name] = [[float(out_df[out_df[0] == i]["Y_test_predicted_proba"]), ind]]
    #print(data)

    rows = []
    cols = []
    for name in list(data.keys()):
        #print(i, name)
        y, x = zip(*data[name])
        rows.append(y)
        for i in x:
            if i not in cols:
                cols.append(i)

    df3 = pd.DataFrame(rows)
    df3.columns = cols
    #df3

    df_init = pd.read_csv(featuredir + "/headerStr_names_init.tab", header=None, sep="\t")
    df_init.columns = ["Name", "Seq"]
    #df_init

    res = pd.concat([df_init, df3], axis=1, sort=False)
    res
    pd.DataFrame.to_csv(res, results + "/out_df_to_rows_TF_{TF_name}_threshold_{threshold}_mono_di.csv".format(TF_name=TF_name, threshold=threshold),  sep='\t', header=True)
    print("Done!")




=======
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
>>>>>>> 51296126fcc189628509f4fc3e053fc75b3d4dbf
