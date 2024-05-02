#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from normalization_PSMSHandler import *
from FDR_LibSearch import *
from mainSearchFunctions import *
import re, numpy as np
import math
import configparser
config = configparser.ConfigParser()
import os
import matplotlib.pyplot as plt
import sys
import time
# In[2]:

start_time = time.time()

# # Read the parameter file using configparser
# params_file = "/Users/spoudel1/Desktop/JUMP_specLib/Program/finalProgram/specLib_firstSearch.params"#sys.argv[1]
params_file = sys.argv[1]

config.read(params_file)

exp_ms2 = config["specLib"]["exp_ms2"] #raw file
specLibFolder = config["specLib"]["specLibFolder"] #library containing folder
ms1_tol = float(config["specLib"]["ms1_tol"]) #precursor ion tolerance
top_ions = int(config["specLib"]["top_ions"]) #total number of top library ions that should be present in spectrum 
ms2_tol = float(config["specLib"]["ms2_tol"]) #fragment ion tolerance

binsize = config["specLib"]["binsize"] #binning width for selecting top ions in the ms2 spectrum 
top_ions_per_bin = int(config["specLib"]["top_ions_per_bin"]) #maximum fragment ions that are selected from raw spectrum for each bin

outputFolder = config["specLib"]["outputFolder"] #search results are stored in this folder

FDR = config["specLib"]["FDR"] #FDR at PSMS level Target Decoy

#tolerace type for dynamic or static tolerance

tolerance_type = config["specLib"]["tolerance_type"]

#Dynamice tolerance file if dynamic tolerance is selected

dyn_tol_file = config["specLib"]["dyn_tol_file"]


#tolerance for ms2 matching = fragment ion tolerance
tol = float(config["specLib"]["tol"])
#example pepxml file to parse the modification information
pepxml = config["specLib"]["pepxml"]

jump_modAA_dict, jump_mod_dict, sta_AA = getDynStatModsInfoPepXml(pepxml)

# exp_ms2 = "/Users/spoudel1/Desktop/JUMP_specLib/Program/ms2/preprocess/FTLD_Batch2_F50.ms2"
specLib = specLibFolder+"/SpectralLibraryTargetDecoy.spLib"
# specLib = "/Users/spoudel1/Desktop/JUMP_specLib/specLib/OneFraction/SpectralLibraryTargetDecoy.spLib"


expMZXML = exp_ms2.split("/")[-1].split(".")[0]


int_sd_dict = {}
if tolerance_type.upper() == "DYNAMIC":
    int_sd_dict = parseDynamicIntensityFile(dyn_tol_file,tol)


# In[5]:


expDF = ms2ToDf_spec(exp_ms2)
libDF = ms2ToDf_spec(specLib)
normalizeIntensity(expDF)
normalizeIntensity(libDF)

logTransformMS2Intensity(expDF)#this is for dynamic intensity tolerance

#library searching results are kept as a dictionary final result
afterCorrQC(expMZXML,expDF,libDF,top_ions,ms2_tol,ms1_tol,top_ions_per_bin,binsize,tolerance_type,int_sd_dict)
