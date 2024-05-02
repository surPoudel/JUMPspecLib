import pandas as pd
from advancedFilter import *
from TargetDecoy import *
from inferenceRT import *
from publicationTables import *
from logFunctions import *
import time
import sys, os, re
import configparser
config = configparser.ConfigParser()

start_time = time.time()

# # Read the parameter file using configparser
# params_file = "/Users/spoudel1/Desktop/JUMP_specLib/Program/finalProgram/specLib_firstSearch.params"#sys.argv[1]
params_file = sys.argv[1]
config.read(params_file)

exp_mzxml = config["specLib"]["exp_mzxml"] #raw file

search_result_file = config["specLib"]["search_result_file"] #search results are stored in this folder

mode = config["specLib"]["mode"]

outputFolder = config["specLib"]["outputFolder"] #folder for fitlering results
#FDR at PSMS level Target Decoy (rt based FDR)
FDR = float(config["specLib"]["initial_FDR"])

### Filter Parameters ###
#User defined Filter method 0 = Rapid (default JDscore), 1 = Rapid (dRT), 10= Advanced filter (deltRT ,deltMZ and Jdscore) 

filter_method = config["specLib"]["filter_method"]

#advanced filtering
binsize = int(config["specLib"]["binsize"])
d_rt = float(config["specLib"]["d_rt"])
d_mz = float(config["specLib"]["d_mz"])
jd = float(config["specLib"]["jdscore"])

#input files for protein level filtering
pitFile = config["specLib"]["pitFile"] #"human_ft_mc2_c57_TMT_K229.pit"
# peptideFile = "idwDecoy_uni_pep.txt"
id_protFile = config["specLib"]["id_protFile"] #"L_ID_mapped_protein_peptides.txt"

###program###

#set current directory as working directory
path = os.getcwd()
os.chdir(path)

##### Protein level filtering required modules #####
#pitrankDict
rankDict = pitFileToRankDict(pitFile)
# id_prot_lookup file to df
idProtDF_noDup = id_prot_lookupFileParse(id_protFile)
selectHighestRankProt(idProtDF_noDup, rankDict)
#peptide and ranked protein dictionary
L_ID_uniqProtDict = dict(zip(idProtDF_noDup.Peptide_ID, idProtDF_noDup.Pit_rank1_protein))
#peptide and unique peptide or not tracker. If unique
#maps to only one protein so value 1 is assigned else np.nan is assigned
uniquePeptideTrack = dict(zip(idProtDF_noDup.Peptide_ID, idProtDF_noDup.uniquePeptideCounter))

######################


#removing log file if previously present
try:
    rmFile("jump_lib_f.log")
except:
    print ("No jump_lib_f.log is present. Writing log file now")
jump_f_outFolder = "sum_"+outputFolder

if mode = 1:
    makedirectory(jump_f_outFolder)

#loads the search result as pandas dataframe
# rank1File = searchResultFolder+"/Library_Search_Rank1.xlsx"

# printDF2Rank1 = pd.read_excel(rank1File)

printDF2Rank1 = mergeSearchResults(search_result_file)



write_log("\n\n****** JUMP -lib filter program *********\n")

write_log("Initiating the JUMP -lib filter program\n")

#if Rapid filteration is selected 
#cutoff description
#absolute delta rt = minutes
#absolute delta mz = ppm
#jd = JDscore 
# d_rt, d_mz, jd = 5, 2, 0.8 

reqd_columns = ['L_ID','Peptide_ID','charge','exp','scan','prec_MZ','Prec_MZ_theoretical','[M+H]+','RT','abs_delRT','JDscore','abs_dPrecMZ','pepLength','Type']

if filter_method != "1":
    if filter_method == "0":
        filter_col = "JDscore"
    if filter_method == "1":
        filter_col = "abs_delRT"

    write_log("\n---- Rapid Filtering method selected ----\n")
    write_log("Input PSMs for filtering = {}\n".format(printDF2Rank1.shape[0]))
    filteredDF_rapid = FDR_Target_Decoy(printDF2Rank1,filter_col)
    filteredScans = publicationTablesRapidFilter(filteredDF_rapid, jump_f_outFolder, exp_mzxml, float(FDR),reqd_columns, L_ID_uniqProtDict, uniquePeptideTrack)
    #copy params file to the new folder
    
else:

    #extraction of RT ms2 based 
    write_log("\n---- Advanced Filtering method selected ----\n")
    #RT_inference_MS2based(exp_mzxml, printDF2Rank1, searchResultFolder)
    
    #RT based filtering and binning
    # psms_firstSearch = binning_FDR_calculation(printDF2Rank1, outputFolderRev, rt_fdr)
    write_log("\nThe binwise report is now being generated ..........\n")

    psms_firstSearch = binning_FDR_calculation(printDF2Rank1, jump_f_outFolder, FDR, d_rt, d_mz, jd, binsize)
    filteredScans = finalPublicationTables(psms_firstSearch, jump_f_outFolder,  exp_mzxml, reqd_columns, L_ID_uniqProtDict, uniquePeptideTrack)

    target = printDF2Rank1.loc[~printDF2Rank1.Peptide_ID.str.contains("Decoy")]
    decoy = printDF2Rank1.loc[printDF2Rank1.Peptide_ID.str.contains("Decoy")]


    write_log("\nGenerating and saving Target Decoy Histograms\n")
    histogramPlot(target, decoy, "JDscore", jump_f_outFolder+"/Target_Decoy_JDscore","target", "decoy")
    histogramPlot(target, decoy, "abs_delRT", jump_f_outFolder+"/Target_Decoy_deltaRT","target", "decoy")
    histogramPlot(target, decoy, "abs_dPrecMZ", jump_f_outFolder+"/Target_Decoy_delMZ","target", "decoy")

#copy params file to the new folder
cmdDir1 = "cp "+params_file+" "+jump_f_outFolder
os.system(cmdDir1)


write_log("The filter program is complete and the associated publication files are generated")
write_log("Time taken for library search and filter --- %s seconds ---{}\n".format(time.time() - start_time))

#move log file to new folder
cmdDir2 = "mv jump_lib_f.log "+jump_f_outFolder
os.system(cmdDir2)