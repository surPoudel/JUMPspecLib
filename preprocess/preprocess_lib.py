#!/usr/bin/env python
# coding: utf-8


import sys
import pandas as pd
import os
import numpy as np
import re
import pyteomics
from pyteomics import mass
from preprocess_functions import *
from job_submission import *
import glob
import configparser
import warnings

from os.path import dirname

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')

config = configparser.ConfigParser()

params_file = os.path.join(os.getcwd(), sys.argv[1])
mzxmls = sys.argv[2]

config.read(params_file)

presearch_path = config["Preprocess"]["presearch_path"] 
pepxml = config["Preprocess"]["pepxml"] 
idtxt = config["Preprocess"]["idtxt"] 
tmtReport = config["Preprocess"]["tmtReport"] 

ion_type_test = config["Preprocess"]["ion_type_test"]
ion_loss_test = config["Preprocess"]["ion_loss_test"] 
resultsDirectory = config["Preprocess"]["resultsDirectory"] 
cluster = config["Preprocess"]["cluster"] 

mkdirCmd = "mkdir "+resultsDirectory
try:
    os.system(mkdirCmd)
except:
    pass


if presearch_path == "0":
    presearch_path = os.getcwd()
    program_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "v2.3/DeisotopeMS1.py")
    preprocess_cmd = "python " + program_path + " {} {}".format(params_file, mzxmls)
    os.system(preprocess_cmd)



# pepxml = "/home/zyuan1/TMT16_SunHuan/Human_FTLD_Rosagrp/1s_1by1/FTLD_Batch2_F1/FTLD_Batch2_F1.1.pepXML"
# idtxt = "/research_jude/rgs01_jude/groups/penggrp/projects/Proteomics_Spectral_Lib/penggrp/step1_generate_library/LIBRARY/SpectralLibrary_TMT16_v0.0.3/HUMAN/jump_f/Reference_human120/sum_HH_tmt16/ID.txt"
# tmtReport = "TMT16"
# resultsDirectory = "Results_15PPM"

# ion_type_test = "a,b,y"
# ion_loss_test = "H2o,NH3"


aminoAcidList = pyteomics.parser.std_amino_acids

ms2_Files_all = glob.glob("{}/*/*.ms2".format(presearch_path))

#remove raw.ms2  files
ms2_Files = []
for ms2file in ms2_Files_all:
    if ".raw.ms2" in ms2file:
        continue
    else:
        ms2_Files.append(ms2file)



dftxt = pd.read_csv(idtxt,delimiter=";",skiprows=return_skiprows(idtxt,";", "Peptide"))
dftxt["spectrum"] = dftxt.apply(createOutfile, df=dftxt, axis=1)
dftxt.drop_duplicates(subset="spectrum", inplace=True, keep="first")
dftxt[["exp","scan","charge"]] = dftxt["spectrum"].str.split(".",expand=True)

dftxt.to_pickle(resultsDirectory+"/dftxt.pkl")


# df_mz = pd.DataFrame(reader)
# protonAdd = mass.calculate_mass(formula='H+')
# df_mz["num"] = df_mz.apply(lambda x: x.params["scan"][0], axis=1)
# df_mz["precursorMz"] = df_mz.apply(lambda x: x.params["precursor m/z"], axis=1)
# df_mz["retentionTime"] = df_mz.apply(lambda x: x.params["RetTime"], axis=1)
# df_mz["CalcPrecNeutralMass"] = df_mz.apply(lambda x: x.params["neutral mass"][0]-protonAdd, axis=1)
# df_mz["charge"] = df_mz.apply(lambda x: x.params["charge"][0], axis=1)




#['num', 'msLevel', 'm/z array', 'intensity array', 'retentionTime', 'precursorMz']




# # Ions that have heavy carbon

# # Define TMT ions for different kinds of TMT




#[exp_mz_list,intensity_list,true_ions_list_parent,matched_int_list_parent,true_ions_ppm_list_parent,matched_theo_ion_list_parent, ionTypes]
#precMZ,monoThPrecursorMass,monoPeptide] #monoisotopic mass and neutral mass of peptide added

# path for the script
source_path_script = sys.argv[0]

#get the path of script
source_path = dirname(source_path_script)

pklfile = resultsDirectory+"/dftxt.pkl"
jobNumbers = []
for in_ms2file in ms2_Files:
    if cluster == "0":
        main(in_ms2file, pklfile, tmtReport, ion_type_test, ion_loss_test, pepxml, resultsDirectory)
    else:
        jobfile = create_job_file(source_path, in_ms2file, pklfile, tmtReport, ion_type_test, ion_loss_test, pepxml, resultsDirectory)
        jobNumber = submit_job(jobfile)
        jobNumbers.append(jobNumber)
        
checkJobStatus(jobNumbers)






'''

# Other parameters
# 0 = disable; 1 = enable; using master node only or entire cluster
#cluster = 0
# LSF used by current cluster; other systems (e.g. SGE & PBS) may be used with minor tweak in source code
#Job_Management_System = LSF
'''


