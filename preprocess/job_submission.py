import os, sys
from preprocess_functions import main

import subprocess
import re
import utils
import time



def create_job_file(source_path, in_ms2file, pkltxt, tmtReport, ion_type_test, ion_loss_test, pepxml, resultsDirectory):
    fileroot = in_ms2file.split("/")[-1].split(".")[0]
    #dta_to_ms2(reqd_dta, mzxml_file)
    log_dir = resultsDirectory+"/log"
    cmd1 = "rm -r "+log_dir
    try:
        os.system(cmd1)
    except:
        print ("No log directory for ", sample)
    os.system("mkdir "+log_dir)
    job_header = "#!/bin/bash\n#BSUB -P preprocess\n#BSUB -J lib_preprocess\n#BSUB -oo {}/{}_log.out\n#BSUB -eo {}/{}_error.err\n#BSUB -n 1\n".format(log_dir, fileroot,log_dir, fileroot)

    preprocess_program = "python {}/main.py {} {} {} {} {} {} {}".format(source_path, in_ms2file, pkltxt, tmtReport, ion_type_test, ion_loss_test, pepxml, resultsDirectory)
   
    jobfile = fileroot+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+preprocess_program)
    return jobfile

def submit_job(jobf):
    cmd = 'bsub -q standard -R "rusage[mem=8000]" < '+jobf
    # os.system(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    (pout, _) = p.communicate()
    pout = pout.decode().strip()
    pout = re.search("(?<=Job <)([0-9]+)", pout)
    jobNumber = pout.group(0)
    return jobNumber

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
