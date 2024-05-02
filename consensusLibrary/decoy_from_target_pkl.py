import pandas as pd
import os, sys, glob
import numpy as np
import pickle
import fileinput
from os.path import dirname
import configparser
from logFunctions import *
from datetime import datetime
import time 
from consensusDecoy import *
from idtxtMs2ModsFunctions import ms2ToDf_spec,mkdir
start_time = time.time()

config = configparser.ConfigParser()

# # Read the parameter file using configparser

params_file = sys.argv[1]
# mzFILE = sys.argv[2]

config.read(params_file)

#Output folder name that contains spectral library files
specLibFolder = config["specLib"]["output_specLibFolder"]
#decoy generation distance
distanceDecoy = float(config["specLib"]["distanceDecoy"])
#1= tmt18_default, 2 = tmt11_default, 3= tmt18_pho, 4 = tmt11_pho, 5= tmt18_ub, 6 = tmt11_ub, 100 = labelfree_default, 1000 = silaclys_default

libtype = float(config["specLib"]["libtype"])

#this dictionaryneeds to be updated with other datatype
libtypeDict = {1:"tmt18_default", 2:"tmt11_default", 3:"tmt18_pho", 4:"tmt11_pho", 5:"tmt18_ub", 6:"tmt11_ub", 100:"labelfree_default", 1000:"silaclys_default"}
libtypename = libtypeDict[int(libtype)]

decoy_gen_method = config["specLib"]["decoy_gen_method"]
targetLib = config["specLib"]["target_library"] #specLibFolder+"/intermediate/jumplib_human_tmt18_default_target.pkl"

# print (decoy_gen_method)
if ((decoy_gen_method != "1") and (decoy_gen_method != "0")):
    print ("Please provide valid decoy generation strategy. Select 0 for precursor swap technique and 1 for mass shift")
    sys.exit(1)

path = os.getcwd()
os.chdir(path)

# make speclib folder
mkdir(specLibFolder)

#make intermediate folder to save intermediate library files to 
# make speclib folder
mkdir(specLibFolder+"/intermediate")


targetDF = pd.read_pickle(targetLib)

#generate the decoy library using newly concatenated library

total_entries = targetDF.shape[0]

appendLibraryFinal_sorted = targetDF.sort_values(['precursorMZ','charge'], ascending=[True, True])
decoy_master_list = [] #decoy_master_dict[key1] = [decoy_scan, precMZ, charge, massNeutral, L_ID, L_peptide, RT,mz];decoy_master_dict[key2] = [decoy_scan_pair, precMZ_pair, charge, massNeutral_pair, L_ID_pair, L_peptide_pair, RT_pair,mz_pair]

if decoy_gen_method == "0": # 0 is precursor swap (default)

    total_blocks = int(total_entries/12000)

    # shuffled = appendLibraryFinal_sorted.sample(frac=1)
    result = np.array_split(appendLibraryFinal_sorted, total_blocks) 
    #decoySpecLibrary_Prec_Swap_New

    exclusion_list = []
    rescue_scans_list = []

    shape = appendLibraryFinal_sorted.shape[0]

    print ("..... Working for {} blocks of target library. Each block has ~12000 scans".format(total_blocks))

    for i,x in enumerate(result):
        print ("..... Working for {} index out of {} blocks of target library".format(i, total_blocks))
        # print ("The section of dataframe columns are {}, and shape of dataframe is {}".format(x.columns, x.shape[0]))
        exclusion_list, rescue_scans_list, decoy_dict = decoySpecLibrary_Prec_Swap_New(x, specLibFolder, distanceDecoy,libtypename, shape, exclusion_list, rescue_scans_list)
        decoy_master_list.append(decoy_dict)



    #check if all scans in rescue_scans_list is in exclusion_list
    rescue_scans_list = list(set(rescue_scans_list).difference(set(exclusion_list)))
    rescue_df_sorted = appendLibraryFinal_sorted[appendLibraryFinal_sorted["scan"].astype("int").isin(rescue_scans_list)]
    rescue_scans_list2 = []

    shape = appendLibraryFinal_sorted.shape[0]
    exclusion_list, rescue_scans_list2,decoy_dict = decoySpecLibrary_Prec_Swap_New(rescue_df_sorted, specLibFolder, distanceDecoy,libtypename, shape, exclusion_list, rescue_scans_list2)
    decoy_master_list.append(decoy_dict)


    decoy_dict = rescue_scan_decoy(appendLibraryFinal_sorted, rescue_scans_list2, distanceDecoy)
    decoy_master_list.append(decoy_dict)
    write_decoy_library(decoy_master_list, specLibFolder, distanceDecoy,libtypename)

else:
    # 1 is precursor mass shift distanceDecoy is the mass shift that is used to generate decoy
    rescue_scans_list2 = list(appendLibraryFinal_sorted.scan.astype("int"))
    decoy_dict = rescue_scan_decoy(appendLibraryFinal_sorted, rescue_scans_list2, distanceDecoy)
    decoy_master_list.append(decoy_dict)
    write_decoy_library(decoy_master_list, specLibFolder, distanceDecoy,libtypename)



# specLib = specLibFolder+"/intermediate/jumplib_human_{}_target.splib".format(libtypename)
#all_prot_ppml = specLibFolder+"/intermediate/id_all_pep.ppml"
specLib = targetLib.split(".")[0]+".splib"
#we need these intermediate files for searching
ppml_intermediate_file = dirname(targetLib)+"/id_all_pep.ppml"

cmd0 = "cp "+ppml_intermediate_file+" "+specLibFolder+"/intermediate"
os.system(cmd0)

# decoyLib = specLibFolder+"/intermediate/jumplib_human_{}_decoy.splib".format(libtypename)
decoyLib_list = glob.glob(specLibFolder+"/intermediate/jumplib_human_*_decoy.splib")
outfilename = specLibFolder+"/jumplib_human_{}.splib".format(libtypename)
filenames = [specLib]+decoyLib_list

#concatenate two spectral library
with open(outfilename, 'w') as fout, fileinput.input(filenames) as fin:
    for line in fin:
        fout.write(line)



cmd2 = "cp "+params_file+" "+specLibFolder
os.system(cmd2)

logFile = "jump_lib_db.log"
cmd3 = "mv "+logFile+" "+specLibFolder
os.system(cmd3)

write_log ("Generating library as the dataframe and store it as a pickle file. This is very important as parsing pickle file can save >1 minute per file during searching.")
libDF = ms2ToDf_spec(outfilename)
libDF.to_pickle(specLibFolder+"/jumplib_human_{}.pkl".format(libtypename))

write_log ("Time taken for generating spectral library database --- %s seconds ---" % (time.time() - start_time))

