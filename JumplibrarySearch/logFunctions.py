import pandas as pd
import os, sys
import multiprocessing as mp
import numpy as np

def write_log(*args):
    with open(args[0], "a") as log_file:
        line = ' '.join([str(a) for a in args[1:]])
        log_file.write(line+'\n')
        print(line)

def rmFile(filename):
    cmd1 = "rm "+filename
    os.system(cmd1)


def makedirectory(folName):
    #create search output directory
    cmdDir = "mkdir "+folName
    
    try:
        os.system(cmdDir)
    except:
        print ("Directory exist")

def fileToDF(file):
    df = pd.read_csv(file, delimiter="\t")
    return df

def parallelize_dataframe(fn, df, args, n_cores=4):
    pool = mp.Pool(n_cores)
    list_of_dict = []
    
    df_split = np.array_split(df, n_cores)
    
    for x in df_split:
        
        ls1 = pool.apply_async(fn, (x,)+args)
        list_of_dict.append(ls1)
    pool.close()
    return list_of_dict

