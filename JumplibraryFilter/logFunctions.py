import pandas as pd
import os, sys

def write_log(*args):
    with open("jump_lib_f.log", "a") as log_file:
        line = ' '.join([str(a) for a in args])
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
