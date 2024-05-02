import pandas as pd
import numpy as np
import sys, os, re
import pickle
import fileinput
from consensusDecoy import *
from DatabaseMergingFunctions import *
from idtxtMs2ModsFunctions import *
from logFunctions import *
from merge_ppml import *
from lowess import genLowessFunction

import configparser
config = configparser.ConfigParser()

import time
import glob

start_time = time.time()


# # Read the parameter file using configparser

params_file = sys.argv[1]
# mzFILE = sys.argv[2]

config.read(params_file)

#unimod modification information
refLib = config["specLib"]["refLib"]
newLib = config["specLib"]["newLib"]
specLibFolder = config["specLib"]["output_specLibFolder"]
distanceDecoy = float(config["specLib"]["distanceDecoy"])

#reference ppml file for all protein
ref_ppmlfile = config["specLib"]["ref_ppmlfile"]

#new ppml file for all protein
new_ppmlfile = config["specLib"]["new_ppmlfile"]

libtype = float(config["specLib"]["libtype"])



#this dictionaryneeds to be updated with other datatype
libtypeDict = {1:"tmt18_default", 2:"tmt11_default", 3:"tmt18_pho", 4:"tmt11_pho", 5:"tmt18_ub", 6:"tmt11_ub", 100:"labelfree_default", 1000:"silaclys_default"}
libtypename = libtypeDict[int(libtype)]

decoy_gen_method = config["specLib"]["decoy_gen_method"]

if (decoy_gen_method != "1") and (decoy_gen_method != "0"):
    print ("Please provide valid decoy generation strategy. Select 0 for precursor swap technique and 1 for mass shift")
    sys.exit(1)

#set current directory as working directory
path = os.getcwd()
os.chdir(path)

mkdir(specLibFolder)
remLog = "rm jump_lib_db.log"

os.system(remLog)



#gets the three dataframes dfR_T, dfL_T,modelRTDF  as a list 
write_log ("\n\n****** PARSING LIBRARIES TO MERGE *********\n")
write_log ("Merging ppml file to generate a consensus ppml file")

gen_merged_ppml(ref_ppmlfile, new_ppmlfile, specLibFolder)



write_log ("\nReading reference library and new library and storing as the dataframes. This will take time\n")

df_all = parseLib_modelDF(refLib, newLib)

#add a function to count peptide occurrence in batches
# print (df_all[0].columns)
# print (df_all[1].columns)
peptideCntDict, batchCntDict = peptidesPerBatches(df_all[0], df_all[1])

#dfR_T = reference library with target only
#dfL_T = new library with target only
# modelRTDF = new dataframe for modeling (only overlapped peptides dataframe with peptide ID, ReferenceRT and NewRT)
dfR_T, dfL_T, modelRTDF = df_all[0], df_all[1], df_all[2] 

notes = list(set(dfR_T.Lib_Description))[0]


#print functions

write_log ("\n\n****** OVERVIEW OF LIBRARY ENTRIES *********\n")

write_log ("Total entries in Reference Library = ",dfR_T.shape[0])
write_log ("Total entries in New Library = ",dfL_T.shape[0])
write_log ("Total entries unique to New Library = ",dfL_T.shape[0]-modelRTDF.shape[0])
write_log ("Total entries unique to Reference Library = ",dfR_T.shape[0] - modelRTDF.shape[0])
write_log ("Total number of overlapped entries = ",modelRTDF.shape[0])

write_log ("\n")
#Keeps track to regenerate same columns dataframe after calibration
columnRequired = dfL_T.columns


#compute minimum and maximum value for the curve generation
minv = np.min(list(dfR_T.RT)+list(dfL_T.RT))-0.1 #this is to avoid the decimal problem
maxv = np.max(list(dfR_T.RT)+list(dfL_T.RT))+0.1 #this is to avoid the decimal problem

write_log ("\n\n****** RT RANGE INFORMATION PRIOR TO CALIBRATION *********\n")

write_log ("Minimum RT value (min) for Reference Library = ",np.min(dfR_T.RT))
write_log ("Maximum RT value (min) for Reference Library = ",np.max(dfR_T.RT))

write_log ("Minimum RT value (min) for New Library = ",np.min(dfL_T.RT))
write_log ("Maximum RT value (min) for New Library = ",np.max(dfL_T.RT))


write_log ("\n\n*************** LOWESS CURVE INFORMATION ******************\n")

write_log ("\nFor Loess curve generation, RT ranging from ", minv, " and ", maxv, " minutes is considered")


write_log ("Generating Loess curve to allign New library RT with Reference RT\n")
f,xnew,ynew = genLowessFunction(modelRTDF, minv, maxv)

#ALl the new RTs are calibrated by applying the fucntion
dfL_T["calibratedRTs"] = dfL_T.apply(lambda x: f(x.RT), axis=1)

#first get the list of all new entries in new database that were absent in Reference database
write_log ("Obtaining the list of all new entries in new database that were absent in Reference database\n")
newEntries = list(set(dfL_T.L_peptide).difference(set(dfR_T.L_peptide)))

#replace the calibrated RTs where isna() in the reference RT column
#this is very important we will not change the RT for the overlapped entries. 
#If you generate new reference RT here, the unique peptides in the old reference need adjustment too.
#All reference entries RT are kept same but new entries are calibrated RTs
df_L_T_new = dfL_T.loc[dfL_T.L_peptide.isin(newEntries)]
# replace RT with calibrated RTs

df_L_T_new["RT"] = df_L_T_new["calibratedRTs"]


#merging of databases
write_log ("\n\n*************** OVERLAPPED ENTRIES DECISION MAKING (Reference or New library) ******************\n")
write_log ("Evaluating overlapped libraries entries. Performing QC based on JScore and making decision on which entry to retain\n")

db_overlaped_merged = QC_newLib_Update(dfR_T, dfL_T, modelRTDF)

#extract reference specific entries
write_log ("Extracting Reference Specific Entries\n")
dfR_unique = dfR_T.loc[~dfR_T.L_peptide.isin(db_overlaped_merged.L_peptide)]

write_log ("Appending Reference Unique Entries with Updated overlaped entries\n")
appendLibrary_R = dfR_unique.append(db_overlaped_merged, sort=False)

write_log ("Sorting the Entries based on Library ID# to ensure same Libray ID# as reference")
appendLibrary_R_sorted = appendLibrary_R.sort_values(by=["L_ID"], ascending=True)

write_log ("Sorting the New Entries based on Library ID# to ensure consistency with future merging")
df_L_T_new_sorted = df_L_T_new.sort_values(by=["L_ID"], ascending=True)

# change the scan number and L_ID for the new_sorted library
new_start_scan = np.max(appendLibrary_R_sorted.scan)+1
new_scan_list = list(range(new_start_scan,new_start_scan+df_L_T_new_sorted.shape[0]))

df_L_T_new_sorted["scan"] = new_scan_list
df_L_T_new_sorted["L_ID"] = df_L_T_new_sorted.scan.apply(lambda x: "p"+str(x).zfill(7))

write_log ("Finally merging all three population together: Reference specific, Updated overlaped and New Library specific Entries")
appendLibraryFinal = appendLibrary_R_sorted.append(df_L_T_new_sorted, sort=False)



# #calibration to be applied only to newRTs ... keeping the reference constant
# dfR_T["calibratedRTs"] = dfR_T["RT"]

# #merging of databases
# write_log ("\n\n*************** OVERLAPPED ENTRIES DECISION MAKING (Reference or New library) ******************\n")
# write_log ("Evaluating overlapped libraries entries.\n")

# write_log ("Calibrating the New Library RT using the Loess curve\n")




# new_peptide_Df = QC_newLib_Update_keep_all_reference(dfR_T, dfL_T, f)


# #loess_newRTs(modelRTDF, dfL_T) #modelRTDF is used for generatting curve #dfL_T is the entire new RT containing dataframe where calibration is applied

# #first get the list of all new entries in new database that were absent in Reference database
# write_log ("Obtaining the list of all new entries in new database that were absent in Reference database\n")
# newEntries = list(set(dfL_T.L_peptide).difference(set(dfR_T.L_peptide)))

# #replace the calibrated RTs where isna() in the reference RT column
# #this is very important we will not change the RT for the overlapped entries. 
# #If you generate new reference RT here, the unique peptides in the old reference need adjustment too.
# #All reference entries RT are kept same but new entries are calibrated RTs
# df_L_T_new = dfL_T.loc[dfL_T.L_peptide.isin(newEntries)]


# write_log ("Sorting the Entries based on Library ID# to ensure same Libray ID# as reference")
# appendLibrary_R_sorted = dfR_T.sort_values(by=["L_ID"], ascending=True)

# write_log ("Sorting the New Entries based on Library ID# to ensure consistency with future merging")
# df_L_T_new_sorted = new_peptide_Df.sort_values(by=["L_ID"], ascending=True)

# write_log ("Finally merging all three population together: Reference specific, Updated overlaped and New Library specific Entries")
# appendLibraryFinal = appendLibrary_R_sorted.append(df_L_T_new_sorted, sort=False)
# print (appendLibraryFinal.columns)

# #output folder is a parameter

# # specLibFolder = "test"
# #peptideCntDict is for observed peptide count across batches
mergeLibrary(appendLibraryFinal,specLibFolder,notes,libtypename, peptideCntDict, batchCntDict)

# #replace all RT with calibrated RTs ### only new peptide RTs are updated else all ref peptide RT is original 
# appendLibraryFinal["RT"] = appendLibraryFinal["calibratedRTs"]

#distanceDecoy is a parameter
# distanceDecoy = 7

write_log ("The Target Updated Libray is now stored. Working now on Decoy library\n")
#idea is to randomly break this huge target library to 
#df = df.sort_values(['id','date'], ascending=[False, True])


total_entries = appendLibraryFinal.shape[0]

appendLibraryFinal_sorted = appendLibraryFinal.sort_values(['precursorMZ','charge'], ascending=[True, True])
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
    #generate the decoy library using newly concatenated library
    decoySpecLibrary(appendLibraryFinal, specLibFolder, distanceDecoy,libtypename)

    # # 1 is precursor mass shift distanceDecoy is the mass shift that is used to generate decoy
    # rescue_scans_list2 = list(appendLibraryFinal_sorted.scan.astype("int"))
    # decoy_dict = rescue_scan_decoy(appendLibraryFinal_sorted, rescue_scans_list2, distanceDecoy)
    # decoy_master_list.append(decoy_dict)
    # write_decoy_library(decoy_master_list, specLibFolder, distanceDecoy,libtypename)




specLib = specLibFolder+"/intermediate/jumplib_human_{}_target.splib".format(libtypename)
# decoyLib = specLibFolder+"/intermediate/jumplib_human_{}_decoy.splib".format(libtypename)
decoyLib_list = glob.glob(specLibFolder+"/intermediate/jumplib_human_{}_decoy.splib".format(libtypename))
outfilename = specLibFolder+"/jumplib_human_{}.splib".format(libtypename)
filenames = [specLib]+decoyLib_list

write_log ("\n\n*************** CONCATENATE TARGET DECOY ******************\n")

write_log ("Concatenating Target Library with Decoy Library\n")
#concatenate two spectral library
with open(outfilename, 'w') as fout, fileinput.input(filenames) as fin:
    for line in fin:
        fout.write(line)

write_log ("The reference library is now updated with new library")
write_log ("Total final entries (Target Only) in updated Reference Library = ", dfR_T.shape[0]+dfL_T.shape[0]-modelRTDF.shape[0])

cmd2 = "cp "+params_file+" "+specLibFolder
os.system(cmd2)

logFile = "jump_lib_db.log"
cmd3 = "mv "+logFile+" "+specLibFolder+"/jump_lib_db_merge.log"
os.system(cmd3)

write_log ("The merging of two databases is complete. The new reference database is located in ",specLibFolder+"/SpectralLibraryTargetDecoy.spLib along with separate Target Library and Decoy Library\n")

write_log ("Merging Libray ID and protein information from both database in new database folder\n")
#this is required for adding protein information 
# extract_merge_L_ID_protTable(refLib, newLib, specLibFolder)


write_log ("\n\n*************** MERGE COMPLETE ******************\n")

write_log ("Generating library as the dataframe and store it as a pickle file. This is very important as parsing pickle file can save >1 minute per file during searching.")

libDF = ms2ToDf_spec(outfilename)
libDF.to_pickle(specLibFolder+"/jumplib_human_{}.pkl".format(libtypename))


write_log ("Time taken for merging two libraries --- %s seconds ---" % (time.time() - start_time))
