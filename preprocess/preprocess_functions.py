import pandas as pd
import pyteomics
from pyteomics import mass,ms2
import numpy as np
import os, sys
import re
from datetime import datetime

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



def precursorIonCompute(peptideMassWithProton, charge):
    precIon = peptideMassWithProton/charge
    return precIon

def calPrec_mz(row):
    peptide = row.plain_peptide
    mods = row.modifications
    pepMass = mass.calculate_mass(peptide)
    charge = int(row.charge)
    protonAdd = mass.calculate_mass(formula='H+')*charge
    modsMassExtract = mods.split(",")
    addmass = 0.0
    for val in modsMassExtract:
        modVal = val.split("_")[2]    
        addmass+=float(modVal)
    
    calcNeutralMass = pepMass+addmass
    peptideMassWithProton = pepMass+addmass + protonAdd
#     print (pepMass, addmass, protonAdd)
#     print (peptideMassWithProton)

    precIon = precursorIonCompute(peptideMassWithProton, charge)
    return pd.Series([precIon,calcNeutralMass])


def ppmCalc(a, b, tol=10):
    a = float(a)
    b = float(b)
  #   massError = abs(a-b)*1e6/a
    massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
    return float(massError)


def correctionFactorCalc_v2(exp_mz_list, intensity_list, theor_fragments_df,tmt, tol=10, ion_types = ["b","y"], ionLoss =["NH3","H2O"]):
    mzIndex = []
    cols = theor_fragments_df.columns
    ionsCols = []
    for ions in ion_types:
        for ionC in cols:
            if ions.lower()+"+" in ionC:
                ionsCols.append(ionC)
            for ionLosses in ionLoss:
                if ions.lower()+"-"+ionLosses.upper() in ionC:
                    ionsCols.append(ionC) 
        

    checkDF = theor_fragments_df[ionsCols]
    newDF = checkDF.replace("---",np.nan)
    arr = newDF.to_numpy().flatten()
    df2 = newDF.copy().stack().reset_index()
    df2.columns = ['oldIndex','ions','m/z']

    ions_type = list(df2.ions)

    #If we want to add tmt ions, just comment the below line out
    theor_ions_list = list(df2["m/z"]) #comment this out
    #theor_ions_list2 = list(arr[~np.isnan(arr)]) #turn this on
        
    #Add TMT ions
    #theor_ions_list=tmt+theor_ions_list2     #turn this on

    true_ions_ppm_list = []
    true_ion_int_list = []
    true_ions_list = []
    true_theretical_ions_list = []
    true_theo_ion_type = []

    for i,exp_ion in enumerate(exp_mz_list):
        for j,theor_ion in enumerate(theor_ions_list):
            ppm = ppmCalc(theor_ion,exp_ion, tol=tol)
            
            if (ppm <= float(tol)) and (ppm >= (-1*float(tol))) :
                true_ions_ppm_list.append(ppm)
                true_ions_list.append(exp_ion)
                true_ion_int_list.append(intensity_list[i])
                true_theretical_ions_list.append(theor_ion)
                true_theo_ion_type.append(ions_type[j])
                mzIndex.append(i)
    #this is added for the QC part of the library
    matchRatioIons = str(len(true_ions_list))+"/"+str(len(theor_ions_list))
    matchPercentage = len(true_ions_list)/len(theor_ions_list)*100

    return true_ions_list, true_ion_int_list,true_ions_ppm_list,true_theretical_ions_list,mzIndex, true_theo_ion_type, matchRatioIons,matchPercentage


#this function is different in consensus library as we cannot sum modifications there
#as the PTMs are used to extract unimod information
def spectrumToDict(spectrum):
    dict1 = {}
    spectrumCommaSplit = spectrum.split(",")
    for x in spectrumCommaSplit:
        y=x.split("_")
        if y[0] not in dict1.keys():
            dict1[y[0]] = [float(y[2])]
        else:
            dict1[y[0]].append(float(y[2]))
    dict2 = {}

    for key in dict1.keys():
        value = np.sum(list(dict1[key]))
        # value = np.sum(list(set(dict1[key])))
        dict2[key] = value
    return dict2


#This fucntion helps tp accurately idenfy phospho position so that the maximum and minium position can be found to calculate phospho loss
def getPhoPosition(massPosDict):
    phoPos = []
    for key, value in massPosDict.items():
        if "79.96" in str(value):
            phoPos.append(int(key))
    return sorted(phoPos)




def checkAA(pepSeq, check_list):
    update_val = "False"
    for aa in list(pepSeq):
        if aa in check_list:
            update_val = "True"
    return update_val




def ionSeriesIonLossSpeRes(peptide,massPosDict, maxcharge=1,useMod ="Yes"):
    #first check if the mods has phosphorylation if yes get the maximum and minimum position
    phoPos = getPhoPosition(massPosDict)
    minPos = 1    #assigns minPosition which is later updated
    maxPos = len(peptide) #assigns maximum position which is later updated
    if len(phoPos) >= 1:
        minPos = np.min(phoPos) #computes minimum phospho position so tha b and y ions can be computed according for phospho loss
        maxPos = np.max(phoPos) #computes maximum phospho position so that b and y ions can be computed according ly

    if useMod != "Yes": #this checks whether the modificaiton is searched or not, generally this is always Yes for localization 
        massPosDict = {0:0.0} #if modificaiton is no than the dictionary has no information
    h2o = mass.calculate_mass(formula='H2O') #h2o mono mass
    co = mass.calculate_mass(formula='CO') #co mono mass
    nh3 = mass.calculate_mass(formula='NH3') #nh3 mono mass
    xmassMore = co-mass.calculate_mass(formula='H2') #hydrogen mono mass
    proton = mass.calculate_mass(formula='H+') #proron mono mass
    hydrogenAtom = mass.calculate_mass(formula='H') #hydrogen atom mono mass
    hpo3 = mass.calculate_mass(formula='HPO3') #computes hpo3 mono mass
    h3po4 = mass.calculate_mass(formula='H3PO4') #computes h3po4 mono mass

    all_ions_dict = {}#iniitates a dictionary for theoretical ions which is later conveted to dataframe
#possible ion loss or no ion loss
    ionloss = {"":0,"-H2O":h2o,"-HPO3":hpo3,"-H3PO4":h3po4,"-NH3":nh3}
    
#this section computes a,b,c ions    
    for i in range(1, len(peptide)+1):
        addModMass = 0.0
        for massKey in massPosDict.keys():
            if int(massKey) <= i:
                addModMass += float(massPosDict[massKey])
                
        valAddKey(all_ions_dict,"Seq",peptide[i-1])
#         print (peptide[0:i-1])
        for losses in ionloss.keys():
            if losses == "-H2O":
#water loss if the amino acid sequence is STED
                fate = checkAA(peptide[0:i],["S","T","E","D"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i < len(peptide)):
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
            if losses == "-NH3":
#ammonia loss if the amono acid sequence is RKQN
                fate = checkAA(peptide[0:i],["R","K","Q","N"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i < len(peptide)):
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                                            
            if (losses == "-H3PO4") or (losses == "-HPO3"):
                for charge in range(1, maxcharge):
                    if "79.96" in str(massPosDict.values()):
#                         fate = checkAA(peptide[0:i],["S","T","Y"])
#                         if fate == "True":
#phosoho loss if massPosDict have phopsho loss in it. If yes based on minPos and maxPos ions are calculated
                        if (i >= minPos) and (i < len(peptide)):
                            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
                        else:
                            valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
                            valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                                            
            if losses == "":
                for charge in range(1, maxcharge):
                    if i < len(peptide):
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o-co-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"c(-1)"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge-hydrogenAtom)
                        valAddKey(all_ions_dict,"c1"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge+hydrogenAtom)
                        valAddKey(all_ions_dict,"c2"+losses+"+"+str(charge),(mass.calculate_mass(peptide[:i])+(charge*proton)-h2o+nh3-ionloss[losses]+addModMass)/charge+(2*hydrogenAtom))

                    else:
                        valAddKey(all_ions_dict,"b"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"a"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c(-1)"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c1"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"c2"+losses+"+"+str(charge),"---")

#this section computes x,y,z ions
    for i in range(0, len(peptide)):
        
        addModMass = 0.0
        for massKey in massPosDict.keys():
            if int(massKey) > i:
                addModMass += float(massPosDict[massKey])
        for losses in ionloss.keys():    
            if losses == "-H2O":
                fate = checkAA(peptide[i:len(peptide)],["S","T","E","D"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i >=1):
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
            
            if losses == "-NH3":
                fate = checkAA(peptide[i:len(peptide)],["R","K","Q","N"])
                for charge in range(1, maxcharge):
                    if (fate == "True") and (i >=1):
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                    else:
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
            
            
            if (losses == "-H3PO4") or (losses == "-HPO3"):
#                 fate = checkAA(peptide[i:len(peptide)],["S","T","Y"])
                for charge in range(1, maxcharge):
                    if "79.96" in str(massPosDict.values()):
                        if (i < maxPos) and (i >=1):
                            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                        else:
                            valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
            if losses == "":
                for charge in range(1, maxcharge):
                    if (i >=1):
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-ionloss[losses]+addModMass)/charge)
                        valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)+xmassMore-ionloss[losses]+addModMass)/charge)

                        valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+hydrogenAtom)
                        valAddKey(all_ions_dict,"z1"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(2*hydrogenAtom))
                        valAddKey(all_ions_dict,"z2"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(3*hydrogenAtom))
                        valAddKey(all_ions_dict,"z3"+losses+"+"+str(charge),(mass.calculate_mass(peptide[i:])+(charge*proton)-nh3-ionloss[losses]+addModMass)/charge+(4*hydrogenAtom))
                    else:
                        valAddKey(all_ions_dict,"y"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"x"+losses+"+"+str(charge),"---")

                        valAddKey(all_ions_dict,"z"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"z1"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"z2"+losses+"+"+str(charge),"---")
                        valAddKey(all_ions_dict,"z3"+losses+"+"+str(charge),"---")
                 
                                

    df = pd.DataFrame(all_ions_dict)
    return df


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


def each_scan_matching(ms2File, df_dyn, tmt, ion_types, ionLoss, jump_mod_dict, sta_AA):
    run = os.path.basename(ms2File).split(".")[0]
    dftxt_1 = df_dyn.loc[df_dyn.exp == run]
    dftxt_1[["plain_peptide","modifications"]] = dftxt_1.apply(computeModifications, jump_mod_dict=jump_mod_dict,sta_AA=sta_AA, peptide = "Peptide", axis=1)
    dftxt_1[["PrecursorIonMass","CalcPrecNeutralMass"]] = dftxt_1.apply(calPrec_mz, axis=1)

    print ("total PSMs to be analyzed = {} for {} ms2 file".format(dftxt_1.shape[0], run))
    
    phage3plot, phage3plotPrecursoMass={},{}
    
    reader = ms2.IndexedMS2(ms2File)  #reading ms2file file using pyteomics
    df_mz = pd.DataFrame(reader)  #dataframe of the ms2file reader 
    #print (dfMz.columns)
#     protonAdd = mass.calculate_mass(formula='H+')
    df_mz["num"] = df_mz.apply(lambda x: x.params["scan"][0], axis=1)
    df_mz["msLevel"] = 2
    df_mz["precursorMz"] = df_mz.apply(lambda x: x.params["precursor m/z"], axis=1)
    df_mz["retentionTime"] = df_mz.apply(lambda x: x.params["RetTime"], axis=1)
#     df_mz["CalcPrecNeutralMass"] = df_mz.apply(lambda x: x.params["neutral mass"][0]-protonAdd, axis=1)
#     df_mz["charge"] = df_mz.apply(lambda x: x.params["charge"][0], axis=1)
    mz_cols = list(dftxt_1.columns)
    
    np_arr = dftxt_1.to_numpy()
    n=0
    for row in np_arr:
        phaseDictMod = {}
        scan = str(row[mz_cols.index("scan")])
        peptide = row[mz_cols.index("plain_peptide")]
        spectrum = row[mz_cols.index("spectrum")]
        maxCharge = int(row[mz_cols.index("charge")])+1
        mods = row[mz_cols.index("modifications")]
        xcorr = row[mz_cols.index("XCorr")]
        protein = row[mz_cols.index("Protein")]
        peptideJump = row[mz_cols.index("Peptide")]
        mz = df_mz.loc[df_mz.num == scan]['m/z array']
        newDF = df_mz.loc[df_mz["num"] == scan]
        #precMZ = list(newDF["precursorMz"])[0][0]["precursorMz"]
        precMZ = row[mz_cols.index("PrecursorIonMass")]
#         precIntensity = list(newDF["precursorMz"])[0][0]["precursorIntensity"]
        
        intensity = df_mz.loc[df_mz.num == scan]['intensity array']
        exp_mz_list = list(mz.tolist()[0])
        intensity_list = list(intensity.tolist()[0])

        monoThPrecursorMass = row[mz_cols.index("PrecursorIonMass")] #required for consensus library 
        monoPeptide = row[mz_cols.index("CalcPrecNeutralMass")]  #Monoisotopic mass of Peptide with tags
  
        phaseDictMod = spectrumToDict(mods)
        
#         print (peptide,maxCharge,phaseDictMod, scan, peptideJump)
        
        

        if maxCharge == 1:
            checkdf_MOD = ionSeriesIonLossSpeRes(peptide,maxcharge=2,massPosDict=phaseDictMod,useMod ="Yes")
        else:
            checkdf_MOD = ionSeriesIonLossSpeRes(peptide,maxcharge=int(maxCharge),massPosDict=phaseDictMod,useMod ="Yes")
        
        true_ions_list_parent, matched_int_list_parent, true_ions_ppm_list_parent, matched_theo_ion_list_parent,mzIndex, ionTypes, matchRatioIons,matchpercentage = correctionFactorCalc_v2(exp_mz_list,intensity_list, checkdf_MOD, tmt, tol=15, ion_types = ion_types, ionLoss = ionLoss)     
#         raw_mz_match = extractRawMzUsingMzindex(mzIndex, exp_mz_list)
        matched_len = len(true_ions_list_parent)
    
        valuePlotParent = [exp_mz_list,intensity_list,true_ions_list_parent,matched_int_list_parent,true_ions_ppm_list_parent,matched_theo_ion_list_parent, ionTypes, 
                           [precMZ]*matched_len,[monoThPrecursorMass]*matched_len,[monoPeptide]*matched_len, [peptideJump]*matched_len]
        phage3plot[spectrum] = valuePlotParent
#         phage3plotPrecursoMass[spectrum] = [precMZ,monoThPrecursorMass,monoPeptide] #monoisotopic mass and neutral mass of peptide added
        n+=1
        
        if n%1000 == 0:
            print ("Total PSMs completed = {}".format(n))
   
    return run, phage3plot 




def createDf_ms2(dict_of_list):
    matched_ions = []
    matched_intensity = []
    matched_ppm = []
    theoreticalIon = []
    raw_ions = []
    ion_type = []
    spectrum = []
    peptide = []
    #precursor data
    prec_mz = []
    prec_theor_mass = []
    prec_neutral_mass = []
    
    for ionKeys in dict_of_list.keys():
        matched_ions+=dict_of_list[ionKeys][2]
        matched_intensity+=dict_of_list[ionKeys][3]
        matched_ppm+=dict_of_list[ionKeys][4]
        theoreticalIon+=dict_of_list[ionKeys][5]
        ion_type+=dict_of_list[ionKeys][6]
        spectrum += [ionKeys.split("\t")[0]]*len(dict_of_list[ionKeys][3])
        
        #precursor data
        prec_mz+=dict_of_list[ionKeys][7]
        prec_theor_mass+=dict_of_list[ionKeys][8]
        prec_neutral_mass+=dict_of_list[ionKeys][9]
        peptide+=dict_of_list[ionKeys][10]

    ionDFmatch = pd.DataFrame(columns=["spectrum","preprocessedIons","matched_intensity","MassError","theoreticalIon","ion_type",
                                       "precursorMZ","precursor_TheoMZ","precursor_neutral_mass","Peptide"])
    
    ionDFmatch["spectrum"] = spectrum
    ionDFmatch["preprocessedIons"] = matched_ions
    ionDFmatch["matched_intensity"] = matched_intensity
    ionDFmatch["MassError"] = matched_ppm
    ionDFmatch["theoreticalIon"] = theoreticalIon
    ionDFmatch["ion_type"] = ion_type
    
    #"precursor m/z","precursor theoretical m/z","precursor neutral mass
    ionDFmatch["precursorMZ"] = prec_mz
    ionDFmatch["precursor_TheoMZ"] = prec_theor_mass
    ionDFmatch["precursor_neutral_mass"] = prec_neutral_mass
    ionDFmatch["Peptide"] = peptide
    
#     bins = [0,10000,100000,1000000, 10000000,100000000,1000000000,10000000000]
#     ionDFmatch['intRange'] = pd.cut(ionDFmatch.matched_intensity, bins=bins, include_lowest=True)
#     ionDFmatch['right'] = ionDFmatch['intRange'].map(attrgetter('right'))
#     ionDFmatch["intensity"] = ionDFmatch.apply(rightInterval,axis=1)
        
    return ionDFmatch




def library_input_ms2(df, resultsDirectory):
    print ("Generating .ms2 files\n")
    now = datetime.now()
    #print (prev_mod)
    print("now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index
    
    df_11 = df.groupby('spectrum', as_index=False).agg(list) #this creates the lists of all rows that has same spectrum taking spectrum as index

    df_11[["exp","scan","charge"]] = df_11.spectrum.str.split(".", expand=True)
    df_11.scan = pd.to_numeric(df_11.scan, errors='coerce')
    df_12 = df_11.sort_values(by="scan")
#     print (df_12)
    proton = mass.calculate_mass(formula="H+")

    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"
    exp = df_11.exp[0]
    
    new_ms2_file = exp+".theo.ms2pep" #modified the extension name to ms2pep
    
    ms2dir = resultsDirectory+"/ms2"

    mkdirCmd = "mkdir "+ms2dir
    try:
        os.system(mkdirCmd)
    except:
        pass

    msFile = ms2dir+"/"+new_ms2_file

    with open(msFile,"w") as new_ms2:
        new_ms2.write(header_ms2)
        for row in df_12.itertuples():
#             print (row)

            mz_list = row.theoreticalIon 
            mz_int_list = row.matched_intensity
            scan = row.scan
            charge = int(row.charge)
            peptide = row.Peptide[0] #extract peptide information too; this is done in order to modify the ms2 to ms2pep
           
            precursorMz = float(row.precursor_TheoMZ[0])

            massNeutral = (precursorMz*charge) - ((charge-1)*proton)

            new_ms2.write("S\t"+str(scan)+"\t"+str(scan)+"\t"+str(precursorMz)+"\n")
            new_ms2.write("I\tP\t"+peptide+"\n")
            new_ms2.write("Z\t"+str(charge)+"\t"+str(massNeutral)+"\n")

            for index, mz_mass in enumerate(mz_list):
                if np.isnan(mz_list[index]) == False:
                    new_ms2.write(str(mz_list[index])+"\t"+str(mz_int_list[index])+"\n")

    print ("MS2 file for", exp, " fraction is complete\n") 


def main(in_ms2file, pkltxt, tmtReport, ion_type_test, ion_loss_test, pepxml, resultsDirectory):


        
    ions_test = ion_type_test.split(",")
    neutralIonsTest = ion_loss_test.split(",")


    TMT16 = ['126.1277259','127.1247608','127.1310808','128.1281157','128.1344356','129.1314705','129.1377905',
             '130.1348253','130.1411453','131.1381802','131.1445001','132.1415350','132.1478550','133.1448899',
             '133.1512098','134.1482447']

    TMT11 = ['126.1277259','127.1247608','127.1310808','128.1281157','128.1344356','129.1314705',
             '129.1377905','130.1348253','130.1411453','131.1381802','131.1445001']

    TMT10 = ['126.1277259','127.1247608','127.1310808','128.1281157','128.1344356',
             '129.1314705','129.1377905','130.1348253','130.1411453','131.1381802']

    TMT8 = ['126.1277259','127.1247608','127.1310808','128.1344356',
            '129.1314705','129.1377905','130.1411453','131.1381802']

    TMT6 = ['126.1277259','127.1247608','128.1344356',
            '129.1314705','130.1411453','131.1381802'] 


    if tmtReport == "TMT16":
        tmt = TMT16
    if tmtReport == "TMT11":
        tmt = TMT11
    if tmtReport == "TMT10":
        tmt = TMT10
    if tmtReport == "TMT8":
        tmt = "TMT8"
    if tmtReport == "TMT6":
        tmt = TMT6

    dftxt = pd.read_pickle(pkltxt)
    jump_modAA_dict, jump_mod_dict, sta_AA = getDynStatModsInfoPepXml(pepxml)
    run, phage3plot  = each_scan_matching(in_ms2file, dftxt, tmt, ions_test, neutralIonsTest, jump_mod_dict, sta_AA)
    ms2_in_df = createDf_ms2(phage3plot)
    library_input_ms2(ms2_in_df, resultsDirectory)