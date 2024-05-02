import pandas as pd
from advancedFilter import *
from TargetDecoy import *
# from inferenceRT import *
from publicationTables import *
from logFunctions import *
from mulipleFractionsHandler import *
pd.options.mode.chained_assignment = None  # default='warn'

import time
import sys, os, re, glob
import configparser
config = configparser.ConfigParser()

start_time = time.time()

# # Read the parameter file using configparser
# params_file = "/Users/spoudel1/Desktop/JUMP_specLib/Program/finalProgram/specLib_firstSearch.params"#sys.argv[1]
params_file = sys.argv[1]
config.read(params_file)

config.read(params_file)

search_result_path = config["specLib"]["search_result_path"] #search results are stored in this folder

# mode = config["specLib"]["mode"]

outputFolder = config["specLib"]["outputFolder"] #folder for fitlering results
#FDR at PSMS level Target Decoy (rt based FDR)
jdscore_cutoff = float(config["specLib"]["jdscore_cutoff"])

#user defined FDR for cutoff
FDR = float(config["specLib"]["FDR"])
### Filter Parameters ###
#User defined Filter method 0 = Rapid (default JDscore), 1 = Rapid (dRT), 10= Advanced filter (deltRT ,deltMZ and Jdscore) 

filter_method = config["specLib"]["filter_method"]
#grouping size, for example groupsize is defined by delRT interval in y-axis and JDscore interval in x-axis it could be 1 minute and 0.1 JDscore rtwidth, jdwidth = min, score
rt_width_jdscore_width = config["specLib"]["rt_width_jdscore_width"]

rt_width_jdscore_widthsplit = re.split(",| ",rt_width_jdscore_width)
rt_width = float(rt_width_jdscore_widthsplit[0])
jdscore_width = float(rt_width_jdscore_widthsplit[-1])
delRT_SD = float(config["specLib"]["delRT_SD"])
specLibFolder = config["specLib"]["specLibFolder"]

#take all protein ppml file first
all_prot_ppml = specLibFolder+"/intermediate/id_all_pep.ppml"
uni_prot_ppml = specLibFolder+"/intermediate/id_uni_pep.ppml"
#new columns of ppml file Peptides, Protein Group#, Protein Accession #, Protein Description, GN, Fate, group, Protein_sub_group, subgroup
ppmlDFall = ppmlFileReformat(all_prot_ppml)
ppmlDFuni = ppmlFileReformat(uni_prot_ppml)
###program###

#set current directory as working directory
path = os.getcwd()
os.chdir(path)


#removing log file if previously present
try:
    rmFile("jump_lib_f.log")
except:
    print ("No jump_lib_f.log is present. Writing log file now")
jump_f_outFolder = "sum_"+outputFolder


makedirectory(jump_f_outFolder)

#loads the search result as pandas dataframe
# rank1File = searchResultFolder+"/Library_Search_Rank1.xlsx"

write_log("\n\n****** JUMP -lib filter program *********\n")

write_log("Initiating the JUMP -lib filter program\n")


if filter_method == "1":
    all_searchfiles = glob.glob("{}/*/*.allScores.csv".format(search_result_path))
    if len(all_searchfiles) == 0:
        print ("Your search result do not include the RT extraction process. Please check the search parameter rt_extr_align = 1 if you want to use RT power to boost your PSMs, else, change your filter_method = 0 here. This will use JDscore only to perform standard Target/Decoy filtering")
        sys.exit(1)
else:
    all_searchfiles = glob.glob("{}/*/*.1.csv".format(search_result_path))

printDF2Rank1 = mergeSearchResults(all_searchfiles, filter_col= "JDscore", jdscore_cutoff=jdscore_cutoff)
write_log("Input PSMs for filtering = {}\n".format(printDF2Rank1.shape[0]))


reqd_columns = ['Peptide','Protein','Outfile','measuredMH','calcMH','ppm','JDscore','L_ID','RT','group','subgroup','pepLength','Type','unique']

if filter_method == "0":
    write_log("\n---- Rapid Filtering method selected ----\n")

    
    filteredDF_rapid = FDR_Target_Decoy(printDF2Rank1,"JDscore")
    filteredScans = finalPublicationTables(filteredDF_rapid, jump_f_outFolder, float(FDR),reqd_columns, ppmlDFall, ppmlDFuni)
    #copy params file to the new folder

    target = printDF2Rank1.loc[~printDF2Rank1.Peptide.str.contains("Decoy")]
    decoy = printDF2Rank1.loc[printDF2Rank1.Peptide.str.contains("Decoy")]


    histogramPlot(target, decoy, "JDscore", jump_f_outFolder+"/Target_Decoy_JDscore","target", "decoy")

else:

    write_log("\n---- Advanced Filtering method selected ----\n")
    write_log("\nThe entire matrix is now grouped into multiple matrixes with y-axis rtwidth = {} and x-axis jdscore width = {}\n".format(rt_width, jdscore_width))
    filteredDF_advanced = rt_jdscore_filtering(printDF2Rank1, rt_width, jdscore_width, delRT_SD, FDR)
    
    write_log("\nRT assisted filtering is done and the final matrix is generated\n")

    filteredScans = finalPublicationTables(filteredDF_advanced, jump_f_outFolder, float(FDR),reqd_columns, ppmlDFall, ppmlDFuni)
    #copy params file to the new folder

    target = printDF2Rank1.loc[~printDF2Rank1.Peptide.str.contains("Decoy")]
    decoy = printDF2Rank1.loc[printDF2Rank1.Peptide.str.contains("Decoy")]

    write_log("\nGenerating and saving Target Decoy Histograms\n")

    histogramPlot(target, decoy, "JDscore", jump_f_outFolder+"/Target_Decoy_JDscore","target", "decoy")

    
#copy params file to the new folder
cmdDir1 = "cp "+params_file+" "+jump_f_outFolder
os.system(cmdDir1)


write_log("The filter program is complete and the associated publication files are generated")
write_log("Time taken for library search and filter --- %s seconds ---{}\n".format(time.time() - start_time))

#move log file to new folder
cmdDir2 = "mv jump_lib_f.log "+jump_f_outFolder
os.system(cmdDir2)
