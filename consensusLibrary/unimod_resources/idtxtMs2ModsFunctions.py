import pyteomics as pyteomics
import pandas as pd
import os, sys, glob, re
from datetime import datetime
import collections
import numpy as np
from logFunctions import *

#input for this function is id_all_pep.txt
#Peptides   Protein Group#  Protein Accession # Protein Description GN
def peptide_protein_map_library(allPepTxt, jump_mod_dict, sta_AA, specLibFolder):
    allPepTxtdf = pd.read_csv(allPepTxt, delimiter="\t", skiprows=return_rows_nullProgrp(allPepTxt, "\t"))
    
    #For modification or PTMs information on the consensus library
    allPepTxtdf[["plain_peptide","modifications"]] = allPepTxtdf.apply(computeModifications, jump_mod_dict=jump_mod_dict,sta_AA=sta_AA,axis=1)
    allPepTxtdf["massPosDict"] = allPepTxtdf.modifications.apply(lambda x: spectrumToDict(x))
    allPepTxtdf["PeptideSeqWithRealDelMass"] = allPepTxtdf.apply(lambda x: addModValInPepSeq(x.plain_peptide,x.massPosDict), axis=1)
    


    #counts the frequency of peptides in the all table
    allPepTxtdfFreq=allPepTxtdf['PeptideSeqWithRealDelMass'].value_counts().reset_index()
    #Decides whether a peptide is unique or not
    allPepTxtdfFreq['Fate'] = np.where(allPepTxtdfFreq.PeptideSeqWithRealDelMass == 1,'Unique', 'Shared')
    #make a unique Peptide identifier dictionary
    uniqPepDict = dict(zip(allPepTxtdfFreq["index"],allPepTxtdfFreq.Fate))
    allPepTxtdf['Fate'] = allPepTxtdf.PeptideSeqWithRealDelMass.map(uniqPepDict)

    reqdCols = ["Peptides","PeptideSeqWithRealDelMass","Protein Group#","Protein Accession #","Protein Description","GN","Fate"]
    df = allPepTxtdf.copy()[reqdCols]
    df[["Protein_grp","Protein_sub_group"]] = df["Protein Group#"].str.split(".",expand=True)
    df["Protein_sub_group"]=df.Protein_sub_group.astype("int")

    df["Protein_grp_num"] = df["Protein_grp"].str.extract('(\d+)').astype(int)


    #since JUMP -f cannot get all the genes of Tremble database (may be a bug). this part rescues GN of the same group and adds if missing
    dfGN = df.dropna(subset=["GN"])
    grpToGeneDict = dict(zip(dfGN.Protein_grp_num, dfGN.GN))
    df["GN"] = df.Protein_grp_num.map(grpToGeneDict)



    
    dfNR = df.loc[df.groupby(['PeptideSeqWithRealDelMass','Protein_grp']).Protein_sub_group.idxmin()]
    #Generate groupby PeptideSeqWithRealDelMass and protein group 
    

    # dfNR["Protein_grp_num"]=dfNR.Protein_grp_num.astype("int")
    dfNR2 = dfNR.loc[dfNR.groupby('PeptideSeqWithRealDelMass').Protein_grp_num.idxmin()]
    dfNR2["Peptide_ProtGrpKey"] = dfNR2.PeptideSeqWithRealDelMass+"_"+dfNR.Protein_grp

    #make dictionary of rank1 protein with Protein Group

    rankDict = dict(zip(dfNR2.Peptide_ProtGrpKey, dfNR2["Protein Accession #"]))
    # #Decides which protein to keep in unique protein list
    # dfNR2['Representative_Protein'] = np.where(dfNR2.Fate == 'Unique',dfNR2["Protein Accession #"], dfNR2.Protein_grp.map(rankDict))
    dfNR2['Representative_Protein'] = dfNR2["Peptide_ProtGrpKey"].map(rankDict)
    #Peptides PeptideSeqWithRealDelMass   Protein Group#  Protein Accession # Protein Description GN  Fate    Protein_grp Protein_sub_group   Protein_grp_num Peptide_ProtGrpKey  Representative_Protein
    #total columns
    
    #final peptide to protein accession dictionary
    peptDictUniProt = dict(zip(dfNR2.PeptideSeqWithRealDelMass, dfNR2["Protein Accession #"]))
    
    df.to_csv(specLibFolder+"/intermediate/id_all_pep.ppml", sep="\t", index=None)
    dfNR2.to_csv(specLibFolder+"/intermediate/id_uni_pep.ppml",sep="\t", index=None)



    return peptDictUniProt #returns all, unique peptide report and a dictionary with unique peptide linked to unique protein accession


#pitfile to rank
def pitFileToRankDict(pitFile):
    pitDF = fileToDF(pitFile)
    rankDict = dict(zip(pitDF.ProteinName, pitDF.index))
    return rankDict


def selectOneProtein(df, rankDict):
    rank1_prot_list = []
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
#     cnt=0
    for row in np_arr: 
        peptide = row[mz_cols.index("Peptide")]
        protein = rankDict[peptide] #protein representative
        rank1_prot_list.append(protein)
#         cnt+=1
#         print (cnt)
    df["ProteinAccession"] = rank1_prot_list

def selectOneProteinWithIndex(df, rankDict, Protein_accession="Protein"):
    rank1_prot_list = []
    mz_cols = list(df.columns)
    np_arr = df.to_numpy()
#     cnt=0
    for row in np_arr: 
        protein_list = row[mz_cols.index(Protein_accession)]
        peptide = row[mz_cols.index("Peptide")]
        if peptide == "Decoy":
            protein_list = ["Decoy"]

        uniqueProt = list(set(protein_list))
        rank_list = []
        for prot in uniqueProt:
            if prot in rankDict.keys():
                rank_list.append(rankDict[prot])
            else:
                rank_list.append(1000000)
                
        ind = np.argsort([-i for i in rank_list])
        #highest rank is ind[0]
        rank1_prot_list.append(protein_list[ind[0]])
#         cnt+=1
#         print (cnt)
    df["ProteinAccession"] = rank1_prot_list
    


#This function previously take only one value if two modificaitons were seen in same position but now it takes two values
def mixedDictPtmPosDelM(massPosDict1, unimod_mod_infor,final_mods_unimod_Details):
    mixed_dict_ptm_pos_delM = {}
    final_DelM = {}
    for posKey, delMass in massPosDict1.items():
        delMass_split = delMass.split("+")#multiple mods on same site is separated by +
        for x in range(0, len(delMass_split)):
            if posKey not in mixed_dict_ptm_pos_delM.keys():
                mixed_dict_ptm_pos_delM[posKey] = [str(delMass_split[x])+"\t"+",".join(final_mods_unimod_Details[float(delMass_split[x])])]
            else:
                mixed_dict_ptm_pos_delM[posKey].append(str(delMass_split[x])+"\t"+",".join(final_mods_unimod_Details[float(delMass_split[x])]))
    for key in mixed_dict_ptm_pos_delM.keys():
        final_DelM[key]= ",".join(mixed_dict_ptm_pos_delM[key])
    return final_DelM


def addModValInPepSeq(plain_peptide,massPosDict): #massPosDict example is massPosList[0]
    modified_peptide = []
    for index, aa in enumerate(plain_peptide):
        pos = str(index+1)
        if pos in massPosDict.keys():
            massShift = massPosDict[pos]
            aa2 = aa+"("+str(massShift)+")"
            modified_peptide.append(aa2)
        else:
            modified_peptide.append(aa)
    return ("".join(modified_peptide))


def modsForReport(mods, peptide):
    newMods = []
    modSplit = mods.split(",")
    for x in modSplit:
        xSplit = x.split("_")
        pos = xSplit[0]
        modKind = xSplit[1]
        aa = peptide[int(pos)-1]
        if modKind == "S":
            modType = "St"
        else:
            modType = "Dy"
        finalMods = str(pos)+"_"+aa+"_"+str(xSplit[2])+"_"+modType
        newMods.append(finalMods)
    modsPrint = ",".join(newMods)
    return modsPrint


#this function stores the dictionary for all ms2 files in list with spectrum information and precursor mz as key
# and ms2 mz and intensity array as values (list)(each as separate lists)

def ms2fileToDict(ms2_list):
    write_log ("  All ms2 are now stored as dictionary")
    psmsDict = {}
    precDict = {}
    for ms2 in ms2_list:
        write_log ("  Working with {}\n".format(ms2))

        cnt = 0 #keeps account of number of scans in the ms2 file
        with open(ms2,"r") as g:
            sample = ms2.split("/")[-1].split(".")[0] #makes sure that there is no "." inside filename
            for line in g:
                if "S\t" in line: #a scan is found
        #             print (line.strip())
                    mz_list = []
                    int_list = []
                    scan = line.strip().split("\t")[1] # scan is recorded
                    precmzlong = float(line.strip().split("\t")[-1])
                    precmz = str(np.round(precmzlong,7)) #precmz found ##round to 7 decimals

                    cnt+=1 #updates scan number with new Scan found

                if "Z\t" in line:
                    charge = line.strip().split("\t")[1]
                    neutralMasslong = float(line.strip().split("\t")[-1])
                    neutralMass = str(np.round(neutralMasslong,7)) #M+H+ ion #round to 7 decimals
                    key = sample+"."+str(scan)+"."+str(charge)
                    psmsDict[key] = [mz_list,int_list]
                    precDict[key] = precmz
        #             print (key)


                if re.match(r"^\d+.*$",line):
        #             print (line.strip())
                    mz_int = line.strip().split("\t")
                    mz = float(mz_int[0])
                    intensity = mz_int[1]
                    mz_list.append(np.round(mz,7)) #round to 7 decimals
                    int_list.append(int(float(intensity))) #convert to integer from float 
        #             psmsDict[key] = [mz_list,int_list]
        write_log ("  Total scans present in {} file = {}\n".format(ms2, str(cnt)))
        
        psmsDict[key] = [mz_list,int_list]
        precDict[key] = precmz
        write_log ("  Done for ... {}\n".format(ms2))
    return psmsDict, precDict




def sortDictStrKeyToIntAndBackToStr(dict1):
    dict2 = {int(k):v for k,v in dict1.items()}
    dict3 = collections.OrderedDict(sorted(dict2.items()))
    final_dict = {str(k):v for k,v in dict3.items()}
    return final_dict

#QC is now updated for Jscore only with the cutoff of 30 score
#This is new definition by Dr. Peng based on our QC analysis

def QC_keep_throw_spectrum(row):
    xcorr = row.XCorr[0]
    mz_int_pairs = row.mz_int_pairs
    matches = len(row.mz_int_pairs[0][0])
    
    if float(xcorr < 30):
        result = "Throw"
    else:
        result = "Keep"

    # if matches < 3:
    #     result = "Throw"
    # elif (matches < 4) & (xcorr < 20):
    #     result = "Throw"
    # else:
    #     result = "Keep"
    return pd.Series([result,xcorr])






######Other Basic functions ######

def return_rows_nullProgrp(file, delimiter):
    with open(file, "r") as f:
        skiprows = []
        no=0
        for line in f:
            if delimiter not in line:# head info of the table
                skiprows.append(no)
            elif len(line.split(delimiter)[1])==0:# null protein group
                skiprows.append(no)
            no+=1
    return skiprows


def return_skiprows(file, delimiter, peptide):
    with open(file, "r") as f:
        skiprows = 0
        for line in f:
            if peptide+delimiter in line:
                break
            else:
                skiprows+=1
    return skiprows


def createOutfile(row, df):
    columns = list(df.columns)
    if "Outfile" in columns:
        outfile_split = row.Outfile.split(".")
        run = outfile_split[-5].split("/")[-1]
        scan = int(outfile_split[-4])
        charge = outfile_split[-2]
        spectrum = run+"."+str(scan)+"."+str(charge)
    else:
        run = row["Run#"]
        scan = row["Scan#"]
        charge = row["z"]
        spectrum = run+"."+str(scan)+"."+str(charge)
        #   print (spectrum)
    return spectrum



#This fucntion computes modifications with all the static and dynamic represenatation
def computeModifications(row, jump_mod_dict, sta_AA, peptide = "Peptides"):
    mod_peptide = row[peptide] #Peptide sequence predicted by the software
    peptideNoFlank = mod_peptide.split(".")[1] #Flanking peptide
    pepNoFlankList = list(peptideNoFlank) #flanking peptide is converted to the list, so the symbols are also member of list along with amino acids

    plain_peptide_list =[] #iniiates plain peptide formation
    for aa in pepNoFlankList:
        if aa in pyteomics.parser.std_amino_acids: #this looks at the standard amino acids in pyteomics and if the value in list is not amino acid for example * or other symbols, then it will discard those and just adds amino acdis
            plain_peptide_list.append(aa) 
    plain_peptide = "".join(plain_peptide_list) #creates a new string or amino acids
    
    dynAA_count = 0 #counts the number of dynamic modifcaiton in the peptide sequence
    for key in jump_mod_dict.keys(): 
        if key in mod_peptide:
            dynAA_count +=1 #updates dynamic modification if it finds the symbol in modified peptide

    mods = [] #this initiates the peptide modificaion
    if "n" in peptideNoFlank: #looks for dynamic n-term
        dynAA_count-=1 #this updates dynAA_count to 1 value less if n-term is tehre because n term is already dealt here so we can reduce that by 1
        mods.append("1_V_"+str(jump_mod_dict["n"])+"_n") #comet modificaiotn type generation
        peptideNoFlank = peptideNoFlank[1:] #if n-term is found this will now remove n-term from plain peptide for consideration because n-term is already considered
    #find index of other dynamic modification symbols

    mod_symbol_index = [] #initiates a list with index of modified symbols
    for key in jump_mod_dict.keys(): 
        for aa_index, aa in enumerate(peptideNoFlank): #looks for index and amino acid in pepptide Nofranl (it has symbols)
            if aa == key: #if symbol is found
                mod_symbol_index.append(aa_index) #the index is appended
    no_remaining_dyn_mods = len(mod_symbol_index) #lenght of total dynamic modificaiton except n-term

    mod_symbol_index.sort()  #sorts the index in ascending order so that we get correct modifications sites. This is very important when we have more than one dynamic modiification
    # print ("mod symbol _index =", mod_symbol_index)
    iterate = 1
    for mod_no in range(0,no_remaining_dyn_mods):
        aa = peptideNoFlank[mod_symbol_index[mod_no]-1]
        position = mod_symbol_index[mod_no]-iterate+1 #position is always greater than index
        mod_mass = jump_mod_dict[peptideNoFlank[mod_symbol_index[mod_no]]]

        iterate+=1
        new_mod = str(position)+"_V_"+str(mod_mass)
        mods.append(new_mod)
    #static modification

    for sta_key in sta_AA.keys():
        if sta_key == "n":
            sta_mod = "1_S_"+str(sta_AA["n"])+"_n"
            mods.append(sta_mod)
        else:
            plain_pep_list = list(plain_peptide)
            for index, aa in enumerate(plain_pep_list):
                if aa == sta_key:
                    sta_mod = str(index+1)+"_S_"+str(sta_AA[aa])
                    mods.append(sta_mod)
    modifications = ",".join(mods)              
    return pd.Series([plain_peptide,modifications])





#this function is different than validator fucntion as the masses are as list here and not added
def spectrumToDict(spectrum):
    dict1 = {}
    spectrumCommaSplit = spectrum.split(",")
    for x in spectrumCommaSplit:
        y=x.split("_")
        if y[0] not in dict1.keys():
            dict1[y[0]] = [str(y[2])]
        else:
            dict1[y[0]].append(str(y[2]))
    dict2 = {}

    for key in dict1.keys():
        value = "+".join(list(dict1[key])) #for fully tryptic peptides
        # value = "+".join(list(set(dict1[key])))
        dict2[key] = value
    return dict2


#utilizes dynamic modification or static modification dictionaries extracted from pep.xml file
##importing the pickle file for unimod modifications = unimod_mod_infor
def unimodModsDict(jump_modAA_dict,sta_AA,unimod_mod_infor):
    #new dictionary that updates all dynamic and static delta mass as key and unimod description as values
    final_mods_unimod_Details = {}

    #checks search dynamic delta mass in unimod dictionary
    for key,value in jump_modAA_dict.items():
        new_key = round(key,5)
        for key_uni in unimod_mod_infor.keys():
            if new_key == round(key_uni,5):
                final_mods_unimod_Details[key] = (unimod_mod_infor[key_uni])

                
    #checks search static delta mass in unimod dictionary
    for key,value in sta_AA.items():
        new_value = round(value,5)
        for key_uni in unimod_mod_infor.keys():
            if new_value == round(key_uni,5):
                final_mods_unimod_Details[value] = unimod_mod_infor[key_uni]
    
    return final_mods_unimod_Details


#this function extracts dynamic and static modifications using pepxml files and stores them as the dictionary
def getDynStatModsInfoPepXml(pepxml):
    f = open(pepxml,"r") #read the file
    line = f.readline()
    var_AA_mass = {} 
    var_AA_symbol = {} #symbol mass dictionary
    stat_AA_mass = {}
    while "<spectrum_query" not in line: #end reading the file if this is sen
    # look for aminocaid modification keyword so that the modification infomration can be parsed
        if "<aminoacid_modification" in line: #found modification information
            if "symbol=" in line.strip(): #only dynamic modification has symbols as static are fixed

                #check the patter
                pattern = '<aminoacid_modification aminoacid="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)" symbol="(.+?)"/>'

                #match for the patter in the line and store them as the tuples there are 4 matches () this is changed to list with list function
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3,4))

                modAA = mods_Description[0] #first element in list is amino acid 
                varMass = float(mods_Description[1]) #second element in list is mass and it is converted to float
                variable = mods_Description[2] #third element in the list is determination of variable of static modification "Y" is variable and "N" is static
                symbol = mods_Description[3] #Symbol. this is used in the dictionary
                valAddKey(var_AA_mass,  varMass, modAA)
    #             valAddKey(var_AA_symbol, symbol, varMass)
                var_AA_symbol[symbol] = varMass #symbol as key and mass as values in the dictionary
                line = f.readline()
            else:
                #this is for static modification so there is no symbol hence we have only three values in the list
                pattern = '<aminoacid_modification aminoacid="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)"/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3))
                modAA = mods_Description[0]
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                stat_AA_mass[modAA] = varMass
    #             valAddKey(stat_AA_mass, modAA, varMass)
                line = f.readline()

        elif "<terminal_modification terminus" in line: #This is for terminal modification such as N-term or C-term
            if "symbol=" in line.strip():
                pattern = '<terminal_modification terminus="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)".+symbol="(.+?)"/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3,4))

                modAA = mods_Description[0].lower()
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                symbol = mods_Description[3]
                valAddKey(var_AA_mass,  varMass, modAA)
    #             valAddKey(var_AA_symbol, symbol, varMass)
                var_AA_symbol[symbol] = varMass
                line = f.readline()
            else:
                pattern = '<terminal_modification terminus="(\w)" massdiff="(\d+\.\d+)" mass="\d+\.\d+" variable="(\w)".+/>'
                mods_Description=list(re.match(pattern, line.strip()).group(1,2,3))
                modAA = mods_Description[0].lower()
                varMass = float(mods_Description[1])
                variable = mods_Description[2]
                stat_AA_mass[modAA] = varMass
    #             valAddKey(stat_AA_mass, modAA, varMass)
                line = f.readline()
        else:
            line = f.readline()
    return var_AA_mass,var_AA_symbol, stat_AA_mass


def valAddKey(dict1, key, val):
    if key not in dict1.keys():
        dict1[key] = [val]
    else:
        dict1[key].append(val)
    return dict1

def mkdir(outputFolder):
    #create search output directory
    cmdDir = "mkdir "+outputFolder
    try:
        os.system(cmdDir)
    except:
        write_log ("Directory exist")

#read library and save it as to save the parsing time for searching for every mzxml file
def ms2ToDf_spec(ms2File):
    g = open(ms2File,"r")
    lines = g.readlines()
    scan_list = []
    charge_list = []
    MH_list = []
    precursorIon_list = []
    L_ID_list = []
    L_peptide_list = []
    L_Protein_list = []
    RT_list = []
    ms2_mz = []
    ms2_int = []

    mz_list = []

    for line in lines:
        if "S\t" in line:
            if len(mz_list) >= 1:
                ms2_mz.append(mz_list)
                ms2_int.append(int_list)
            mz_list = []
            int_list = []
    #         print ("spectrum found")
            temp_line = line.strip().split("\t")
            scan = int(temp_line[1])
            scan_list.append(scan)
            precursorIon = float(temp_line[-1])
            precursorIon_list.append(precursorIon)
        elif "Z\t" in line:
    #         print ("Charge found")
            temp_line = line.strip().split("\t")
            charge = int(temp_line[1])
            charge_list.append(charge)

            MH = float(temp_line[-1])
            MH_list.append(MH)

        elif "L\tID_with_Modification" in line:
    #         print ("peptide ID found")
            temp_line = line.strip().split("\t")
            L_ID = temp_line[2]
            L_peptide = temp_line[3]
            L_ID_list.append(L_ID)
            L_peptide_list.append(L_peptide)
        elif "L\tProtein" in line:
            temp_line = line.strip().split("\t")
            L_Protein = temp_line[-1]
            L_Protein_list.append(L_Protein)
        
        elif "L\tRT" in line:
    #         print ("RT found")
            temp_line = line.strip().split("\t")
            RT = float(temp_line[2])
            RT_list.append(RT)
        elif re.match(r"^\d+.*$",line):
    #         print ("ms2 FOund")
            #temp_line = line.strip().split("\t")
            temp_line = re.split("\t| ",line) #to adapt to Zuo-Fei's ms2 pattern
            mz_list.append(float(temp_line[0]))
            int_list.append(float(temp_line[1]))

    ms2_mz.append(mz_list)
    ms2_int.append(int_list)

    dict1 = {"scan":scan_list,"charge":charge_list,"[M+H]+":MH_list,
                "prec_MZ":precursorIon_list,"L_ID":L_ID_list,"L_peptide":L_peptide_list,
                "L_protein":L_Protein_list,"RT":RT_list,"m/z":ms2_mz,"intensity":ms2_int}

    dict2 = {}
    for key in dict1.keys():
        if len(dict1[key]) != 0:
            dict2[key] = dict1[key]

    ms2Df = pd.DataFrame.from_dict(dict2)
    
    return ms2Df
