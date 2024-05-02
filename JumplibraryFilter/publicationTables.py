import pandas as pd 
import re
import os, sys
import numpy as np
from TargetDecoy import *
from advancedFilter import *

def finalPublicationTables(psms_firstSearch, outputFolder,  FDR, reqd_columns, ppmlDFall, ppmlDFuni):
    #to publish psms level and peptide level information, we may remove some columns 

    #make the file with similar to JUMP ID.txt at least Peptide column, exp etc
    write_log ("A cutoff of {}% FDR is applied for filtering\n".format(FDR))
    psms_firstSearch3 = psms_firstSearch.query("FDR<{}".format(FDR))
    psms_firstSearch3=psms_firstSearch3.drop_duplicates(subset="Outfile", keep="first")
    
    #peptide to unique protein dictionary
    uniqueProtDict = dict(zip(ppmlDFuni.Peptide, ppmlDFuni.Representative_Protein))

    write_log ("IDs (Target + Decoy) file written in {}/IDwDecoy.txt file\n".format(outputFolder))

    #For publication only
    #protein level publication filtering and publication tables
    # psms_firstSearch3.loc[psms_firstSearch3.Protein.isna(),"Protein"] = psms_firstSearch3["L_ID"]

    psms_firstSearch3["Unique_Protein"] = psms_firstSearch3.Peptide.map(uniqueProtDict)
    #replace NA with Unique_Protein with L_ID. This happens with Decoy matches
    psms_firstSearch3["Unique_Protein"]=np.where(psms_firstSearch3.Unique_Protein.isnull(),psms_firstSearch3.L_ID, psms_firstSearch3.Unique_Protein)

    #replace NA with Protein with L_ID. This happens with Decoy matches
    psms_firstSearch3["Protein"]=np.where(psms_firstSearch3.Protein.isnull(),psms_firstSearch3.L_ID, psms_firstSearch3.Protein)
    
    #peptide to protein dictionary
    #this dictionary is used for protein level publications
    uniqProtAllProtDict = dict(zip(psms_firstSearch3.Unique_Protein, psms_firstSearch3.Protein))

    psms_firstSearch3_all = tidy_split(psms_firstSearch3, "Protein", sep=',', keep=False)
    psms_firstSearch3_all["Protein Accession #"] = psms_firstSearch3_all["Protein"]
    psms_firstSearch3_all_ppml = psms_firstSearch3_all.merge(ppmlDFall, how="left", on=["Peptide","Protein Accession #"])

    #here since ppml file has target only, decoy will be null so we replace the decoy with not unique else keep the original value
    psms_firstSearch3_all_ppml["unique"]=np.where(psms_firstSearch3_all_ppml.unique.isnull(),0, psms_firstSearch3_all_ppml.unique)

    decoy = 0
    try:
        decoy = int(psms_firstSearch3.groupby('Type').count()["scan"].Decoy)
    except:
        print ("....No Decoys found")
    target = int(psms_firstSearch3.groupby('Type').count()["scan"].Target)

    #list of scans that are avoided during second search 
    filteredScans = list(psms_firstSearch3.scan)

    write_log ("ID (Target only) file written in {}/IDwDecoy.txt file\n".format(outputFolder))
    write_log ("Total Target IDs = {}".format(target))
    write_log ("Total Decoy IDs = {}".format(decoy))
    write_log ("Final FDR = {}".format(decoy/target*100))
    

    idtarget = psms_firstSearch3_all_ppml[reqd_columns].loc[~psms_firstSearch3_all_ppml[reqd_columns].Peptide.str.contains("Decoy")]
    idtarget.to_csv(outputFolder+"/ID.txt", sep="\t",index=None)

    psms_firstSearch3_all_ppml[reqd_columns].to_csv(outputFolder+"/IDwDecoy.txt", sep="\t",index=None)
    psms_firstSearch3.to_csv(outputFolder+"/psms_qval.txt", sep="\t",index=None)


    peptideFilteredDF, peptideFilteredDFall  = peptideLevelFiltering(psms_firstSearch3_all_ppml)
    
    peptideTableCols = ["Peptides","Protein Group#","Protein Accession #","Protein Description","GN","PSM#","Run#","Scan#","m/z","z","ppm","JDscore","L_ID","RT","abs_dPrecMZ","unique"]
    pubFolder = outputFolder+"/publications"
    makedirectory(pubFolder)

    #reporting target peptides
    peptideFilteredDF[peptideTableCols].to_csv(pubFolder+"/id_uni_pep.txt", sep="\t", index=None)
    peptideFilteredDFall[peptideTableCols].to_csv(pubFolder+"/id_all_pep.txt", sep="\t", index=None)
    

    
    #protein level FDR; need idDecoy_uni_pep.txt to compute FDR
    proteinFilteredDF_uni, proteinFilteredDF_all = proteinLevelFiltering(peptideFilteredDF, uniqProtAllProtDict)


    publishCols = ["Protein Group#","Protein Accession #","Protein Description","GN","PSM#","Total Peptide#", "Unique Peptide#","Peptide of the Highest Score","Run#", "Scan#","m/z","z","ppm","L_ID","RT","JDscore","abs_dPrecMZ"]

    
    proteinFilteredDF_uni[publishCols].to_csv(pubFolder+"/id_uni_prot.txt", sep="\t", index=None)
    proteinFilteredDF_all[publishCols].to_csv(pubFolder+"/id_all_prot.txt", sep="\t", index=None)

    return filteredScans


