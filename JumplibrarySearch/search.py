#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from LibrarySearchFunctions import *
from FDR_LibSearch import *
import re, numpy as np
import math
import configparser
config = configparser.ConfigParser()
import os
import matplotlib.pyplot as plt



# # Read the parameter file using configparser

params_file = "/Users/spoudel1/Desktop/JUMP_specLib/Program/finalProgram/specLib_firstSearch.params"#sys.argv[1]

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


# In[39]:


# exp_ms2 = "/Users/spoudel1/Desktop/JUMP_specLib/Program/ms2/preprocess/FTLD_Batch2_F50.ms2"
specLib = specLibFolder+"/SpectralLibraryTargetDecoy.spLib"
# specLib = "/Users/spoudel1/Desktop/JUMP_specLib/specLib/OneFraction/SpectralLibraryTargetDecoy.spLib"




int_sd_dict = {}
if tolerance_type.upper() == "DYNAMIC":
    int_sd_dict = parseDynamicIntensityFile(dyn_tol_file,tol)


# In[54]:


expDF = ms2ToDf_spec(exp_ms2)
libDF = ms2ToDf_spec(specLib)
normalizeIntensity(expDF)
normalizeIntensity(libDF)

logTransformMS2Intensity(expDF)#this is for dynamic intensity tolerance



# dfDict = {}#collection of matched precursor MZ dataframe
# ms2_tol = 15 #this is parameter for precursor ion matches
mz_cols = list(expDF.columns)
np_arr = expDF.to_numpy()
# top_ion=3
dot_product_results = {}
cnt = 1
if tolerance_type.upper() == "DYNAMIC":
    ms2_tol = max(int_sd_dict.values())

for row in np_arr:
    scan = str(row[mz_cols.index("scan")])
    charge = int(row[mz_cols.index("charge")])
    precMZ = float(row[mz_cols.index("prec_MZ")])
    mz = row[mz_cols.index("m/z")]
    intensity = list(row[mz_cols.index("normalized_intensity")])
    log10_intensity = list(row[mz_cols.index("log10_intensity")]) #this intensity is used for dynamic tolerance
    
    matched_lib_DF = scanPrecursorMatch(precMZ, charge, libDF, ms1_tol)
    
    if matched_lib_DF.shape[0] >=1: #there may not be any matches to the library so we have to check the library match to avoid error
#         print ("scan number is ",scan)
#         dfDict[scan] = matched_lib_DF
        matched_lib_DF_top2 = checkTopLibraryIons(mz, matched_lib_DF, top_ions,ms2_tol)
        if binsize.upper() == "NO":
            spectrumInputDict = {"mz":mz,"intensity":intensity}
            tolInputDict = {"mz":mz,"intensity":log10_intensity}
        else:
            binsize = int(binsize)
            spectrumInputDict = binning_mz_100(mz, intensity, top_ions_per_bin,binsize) #simplified dictionary with 10 top ions for each 100 mz bins with normalized intensity 
            #this is the input for dot product that will be compared with library
        
            #the input dictionary for tolerance dynamic tolerance
            tolInputDict = binning_mz_100(mz, log10_intensity, top_ions_per_bin,binsize) #SD input file has intensity as log10 so we need log10_intensity
        
        mz_cols2 = list(matched_lib_DF_top2.columns)
        np_arr2 = matched_lib_DF_top2.to_numpy()
        for row2 in np_arr2: 
            lib_mz_int_dict = {}
            L_ID = row2[mz_cols2.index("L_ID")]
            RT = row2[mz_cols2.index("RT")]
            L_peptide = row2[mz_cols2.index("L_peptide")]
#             print(L_peptide)
            lib_mz_int_dict["mz"] = row2[mz_cols2.index("m/z")]
            lib_mz_int_dict["intensity"] = list(row2[mz_cols2.index("normalized_intensity")]) #use normalize dictionary for the library input

#             #tolerance_type defines how to choose the tolerance Dynamic will use the int_sd_dict dictionary else static tolerace will be used
            tr_featSpec = trimFeatSpec(spectrumInputDict,tolInputDict,int_sd_dict,tolerance_type,lib_mz_int_dict) #int_sd_dict = tolerance SD defined intensity from the file after mass calibration
            normalized_dp = calcMS2Similarity(tr_featSpec,lib_mz_int_dict,ms2_tol) #fragment ion tolerance
            if int(scan) not in dot_product_results.keys():
                dot_product_results[int(scan)] = [L_ID+";"+str(normalized_dp)+";"+str(RT)+";"+L_peptide]
            else:
                dot_product_results[int(scan)].append(L_ID+";"+str(normalized_dp)+";"+str(RT)+";"+L_peptide)
#         print ("Total scan searched = ", cnt,"\tscan=" ,scan)
#         cnt+=1
#     else:
#         dot_product_results[int(scan)] = ["NA;0.0;0.0;Decoy_fake_NA"]
final_result = {}
for results in dot_product_results.keys():
#     print (len(dot_product_results[results]))
    final_result[results] = ",".join(dot_product_results[results])
    
    
 

#create search output directory
cmdDir = "mkdir "+outputFolder
os.system(cmdDir)
try:
    os.system(cmdDir)
except:
    print ("Directory exist")



expDF["simMS2"] = expDF.scan.map(final_result)
colKeep = ['scan', 'charge', '[M+H]+', 'prec_MZ','simMS2']
expDF_simMS2 = expDF[colKeep]
expDF_simMS2_dropNA_simMS2 = expDF_simMS2.dropna(subset=["simMS2"])
expDF_simMS2_dropNA = expDF_simMS2_dropNA_simMS2.copy(deep=False)
expDF_simMS2_dropNA["simMS2_Ranked"] = expDF_simMS2_dropNA.apply(rankMatchedPSMS, axis=1)
expDF_simMS2_dropNA_split = tidy_split(expDF_simMS2_dropNA, "simMS2_Ranked", sep=',', keep=False)

#we now split simMS2_Ranked

expDF_simMS2_dropNA_split[["L_ID","Library_Match_score (DP)","RT","Peptide_ID","Charge_check","Prec_MZ_theoretical","Rank (PSMS)"]]=expDF_simMS2_dropNA_split.simMS2_Ranked.str.split(";",expand=True) 

#convert string scores to float scores
expDF_simMS2_dropNA_split["LDscore"] = expDF_simMS2_dropNA_split["Library_Match_score (DP)"].astype("float")

printCols = ['scan', 'charge', '[M+H]+', 'prec_MZ', 'L_ID', 'RT', 'Peptide_ID','Rank (PSMS)','Prec_MZ_theoretical', 'LDscore']
printDF = expDF_simMS2_dropNA_split[printCols]
printDF2 = printDF.sort_values(by=["LDscore"], ascending=False)
printDF2.to_excel(outputFolder+"/First_Library_Search_All_Ranks.xlsx",index=None)
#select rank1 for histogram
printDF2Rank1 = printDF2.loc[printDF2["Rank (PSMS)"] == "Rank1"]
printDF2Rank1.to_excel(outputFolder+"/First_Library_Search_Rank1.xlsx", index=None)


printDF2Rank1FDR = FDR_Target_Decoy(printDF2Rank1)

#remove decoy hits and select IDs only less than FDR=parameter
psms_firstSearch = printDF2Rank1FDR.loc[(~printDF2Rank1FDR.Peptide_ID.str.contains("Decoy")) & (printDF2Rank1FDR.FDR <= 1.0)]

psms_firstSearch.to_csv(outputFolder+"/ID.txt", sep="\t",index=None)


# # plot histogram

target = printDF2Rank1FDR.loc[~printDF2Rank1FDR.Peptide_ID.str.contains("Decoy")]
decoy = printDF2Rank1FDR.loc[printDF2Rank1FDR.Peptide_ID.str.contains("Decoy")]


minv = np.min(decoy.LDscore)
maxv = np.max(target.LDscore)
testbinsize = np.linspace(minv,maxv)


plt.rcParams.update({'font.size': 10})
fig,ax = plt.subplots(figsize=(4,2))
plt.yticks(color="black")
# size, scale = 1000, 10

commutes2 = target['LDscore']
commutes2.plot.hist(grid=False, bins=50, rwidth=0.9,
                   color='#F4F6F7',edgecolor='black', linewidth=1.0)
commutes = decoy['LDscore']
commutes.plot.hist(grid=False, bins=50, rwidth=0.9,
                   color='#808B96',edgecolor='black', linewidth=1.0)
# Hide grid lines
plt.grid(False)

plt.title('')
plt.xlabel('MS2 matching score (Dot Product)')
plt.ylabel('Number of PSMs')
# plt.grid(axis='y', alpha=0.75)
plt.legend(["target","decoy"],loc="best")
figurename = outputFolder+"/Target_Decoy_7Da_AllScans_NULL.pdf"
figurename1 = outputFolder+"/Target_Decoy_7Da_AllScans_NULL.png"
fig.savefig(figurename, bbox_inches="tight", dpi=600 )
fig.savefig(figurename1, bbox_inches="tight", dpi=600 )

