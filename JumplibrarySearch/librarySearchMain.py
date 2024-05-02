import os, sys, glob
from os.path import dirname
import argparse
from logFunctions import *
import configparser
import subprocess
import re
import utils
import time
config = configparser.ConfigParser()



def msg(name=None):
    return '''\n\n\npython mainProgram.py params_file\n\n\n '''

parser = argparse.ArgumentParser(description="JUMP spectral library search", prog="mainProgramServer.py",usage=msg())
parser.add_argument("parameterfile", help="Jump-lib-s parameter file")
parser.add_argument("--queue",dest="queue", default="queue=standard")
parser.add_argument("--mem",dest="mem", default="mem=20000")

args = parser.parse_args()
queue = args.queue.split("=")[-1]
mem = args.mem.split("=")[-1]

params_file = args.parameterfile

config.read(params_file)

# mzxml_path = config["specLib"]["mzxml_path"]
ms2_path = config["specLib"]["ms2_path"]
job_type = config["specLib"]["job_type"]

# path for the script
source_path_script = sys.argv[0]
#get the path of script
source_path = dirname(source_path_script)


def create_job_file(params_file,ms2File):
    fileroot = ms2File.split("/")[-1].split(".")[0]
    makedirectory(fileroot)
    path = os.getcwd()

    log_dir = path+"/"+fileroot+"/log"
    os.system("mkdir "+log_dir)
    job_header = "#!/bin/bash\n#BSUB -P JUMP-lib-s\n#BSUB -J spectral library search\n#BSUB -oo "+log_dir+"/log.out\n#BSUB -eo "+log_dir+"/error.err\n#BSUB -n 1\n"
    job_body1 = "python {}/librarySearch.py {} {}".format(source_path,params_file,ms2File)
    jobfile = fileroot+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+job_body1)
    return jobfile


def checkJobStatus(jobNumbers):
    nRemainder = len(jobNumbers)
    while nRemainder > 0:
        time.sleep(10)
        p = subprocess.Popen("bjobs -noheader", shell=True, stdout=subprocess.PIPE)
        (pout, _) = p.communicate()
        pout = pout.decode().strip().split("\n")
        jobRemainder = [p.split()[0] for p in pout]
        nRemainder = len(set(jobNumbers).intersection(jobRemainder))
        nFinished = len(jobNumbers) - nRemainder
        text = "\r  {} job(s) is/are finished".format(nFinished)
        sys.stdout.write(text)
        sys.stdout.flush()
    print ("  {} job(s) is/are finished".format(nFinished))

def submit_job(jobf,queue,mem):
    cmd = 'bsub -q '+queue+' -R "rusage[mem='+mem+']" < '+jobf
    # os.system(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    (pout, _) = p.communicate()
    pout = pout.decode().strip()
    pout = re.search("(?<=Job <)([0-9]+)", pout)
    jobNumber = pout.group(0)
    return jobNumber

ms2_files_all = glob.glob(ms2_path+"/*.ms2")

if len(ms2_files_all) == 0:
    print ("There is no ms2 file in {} folder. Trying subfolder in this path".format(ms2_path))
    ms2_files_all = glob.glob(ms2_path+"/*/*.ms2")
    if len(ms2_files_all) == 0:
        print ("There are no ms2 files in the folder. Come back with ms2 files in the folder\n")
#need to remove .raw.ms2 from the list
ms2_files = []

for file in ms2_files_all:
    if ".raw." not in file:
        ms2_files.append(file)

if job_type == "0":
    for ms2File in ms2_files:
        cmd =  "python {}/librarySearch.py {} {}".format(source_path,params_file,ms2File)
        os.system(cmd)

else:
    

    jobNumbers = []
    for ms2File in ms2_files:
        job_num = submit_job(create_job_file(params_file,ms2File),queue,mem)
        print ("\n\nJob is submitted for spectral library searching for {}\n".format(ms2File))
        jobNumbers.append(job_num)

    checkJobStatus(jobNumbers)





