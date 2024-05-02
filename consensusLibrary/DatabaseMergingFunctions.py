import pandas as pd
import numpy as np
import sys, os, re
from os.path import dirname
from scipy import interpolate
import statsmodels.api as sm
from datetime import datetime
from consensusDecoy import *
from logFunctions import *

def targetLibToDF(ms2File):
    g = open(ms2File,"r")
    lines = g.readlines()
    scan_list = []
    charge_list = []
    MH_list = []
    precursorIon_list = []
    L_ID_list = []
    L_peptide_list = []
    L_prot_list = []
    Mod_Pos_list = [] #append mod position line
    spectrum_list = [] #append spectrum line
    jscore_list = [] #append jscore line
    psm_no_list = [] #append psm no line
    RT_list = []
    libraryNotes_list = []
    batchCntList = []
    pepCntList = []
    ms2_mz = []
    ms2_int = []

    mz_list = []
    max_Jscore_list = []

    for index,line in enumerate(lines):
        if "S\t" in line:
            if len(mz_list) >= 1:
                ms2_mz.append(mz_list)
                ms2_int.append(int_list)
                Mod_Pos_list.append(mod_line_list)
            mod_line_list =[]
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

            MH = float(temp_line[2])
            MH_list.append(MH)

        elif "L\tID_with_Modification" in line:
    #         print ("peptide ID found")
            temp_line = line.strip().split("\t")
            L_ID = temp_line[2]
            L_peptide = temp_line[3]
            L_ID_list.append(L_ID)
            L_peptide_list.append(L_peptide)
        
        elif "L\tProtein\tRepresentative" in line:
    #         print ("peptide ID found")
            temp_line = line.strip().split("\t")
            L_prot = temp_line[3]
            L_prot_list.append(L_prot)

        elif "L\tMod_Pos" in line:
            
            mod_line_list.append(line.strip())
            
        
        elif "L\tSpectrum" in line:
            spectrum_list.append(line.strip())
        
        elif "L\tJScore" in line:
            jscore_list.append(line.strip())
            maxjscore = float(line.strip().split("\t")[2].split(";")[0])
            max_Jscore_list.append(maxjscore)
            
        elif "L\tPSM#" in line:
            psm_no_list.append(line.strip())
        
         
        elif "L\tRT" in line:
    #         print ("RT found")
            temp_line = line.strip().split("\t")
            RT = float(temp_line[2])
            RT_list.append(RT)

        #new_ms2.write("L\tLibraryNotes"+"\t"+libraryNotes+"\n")
        elif "L\tLibraryNotes" in line:
            temp_line = line.strip().split("\t")
            notes = temp_line[2]
            libraryNotes_list.append(notes)

        elif "L\tTotalBatches" in line:
            temp_line = line.strip().split("\t")
            batchCnt = int(float(temp_line[2]))
            batchCntList.append(batchCnt)


        #"L\tPeptideLibraryCount\t1\n"
        elif "L\tPeptideLibraryCount" in line:
            temp_line = line.strip().split("\t")
            peptideCnt = int(float(temp_line[2]))
            pepCntList.append(peptideCnt)

        elif re.match(r"^\d+.*$",line):
    #         print ("ms2 FOund")
            #temp_line = line.strip().split("\t")
            temp_line = re.split("\t| ",line) #to adapt to Zuo-Fei's ms2 pattern
            mz_list.append(float(temp_line[0]))
            int_list.append(float(temp_line[1]))
    
    #print ("Total scans = ", len(scan_list))
    #print ("Total notes = ",len(libraryNotes_list))
    Mod_Pos_list.append(mod_line_list)
    ms2_mz.append(mz_list)
    ms2_int.append(int_list)

    dict1 = {"scan":scan_list,"charge":charge_list,"[M+H]+":MH_list,
                "precursorMZ":precursorIon_list,"L_ID":L_ID_list,"L_peptide":L_peptide_list,"L_protein":L_prot_list,
             "mod_pos_line":Mod_Pos_list,"spectrum_line":spectrum_list,"jscore_line":jscore_list,"max_JScore":max_Jscore_list,
              "psm_no_line":psm_no_list,"RT":RT_list,"Lib_Description":libraryNotes_list,"Total_Batches":batchCntList,"Peptide_observed":pepCntList,"m/z":ms2_mz,"intensity":ms2_int}



    dict2 = {}
    for key in dict1.keys():
        if len(dict1[key]) != 0:
            dict2[key] = dict1[key]

    ms2Df = pd.DataFrame.from_dict(dict2)
    return ms2Df



def parseLib_modelDF(refLib, newLib):
    write_log ("Parsing reference Library and storing as dataframe")


    try:
        dfR = targetLibToDF(refLib)
    except:
        write_log ("Please use only Target Reference library. Name of library is jumplib_human_tmt18_target.splib")
        sys.exit(1)

    write_log ("Parsing new Library and storing as dataframe")

    try:
        dfL = targetLibToDF(newLib)
    except:
        write_log ("Please use only Target New library. Name of library is jumplib_human_tmt18_target.splib")
        sys.exit(1)


    # dfR["Lib_Description"] = "Reference"
    # dfL["Lib_Description"] = "DongGeun/Boer/Danting TMT16 experiments htmt_b301 batch"
    
    dfR_T = dfR.loc[~dfR.L_peptide.str.contains("Decoy")]
    dfL_T = dfL.loc[~dfL.L_peptide.str.contains("Decoy")]
    reqdCols = ["L_peptide","RT"]
    
    
    overlapTargetsDF = dfR_T[reqdCols].merge(dfL_T[reqdCols], how="inner", on="L_peptide")
    renameColsDict = {"RT_x":"ReferenceRT","RT_y":"NewRT"}
    
    modelRTDF = overlapTargetsDF.rename(columns=renameColsDict)
    
    return [dfR_T, dfL_T, modelRTDF]


#function to add peptide occurrence in multiple batches during merging the library

def peptidesPerBatches(df1, df2):
    dfMerge = df1.merge(df2, how="outer", on="L_peptide")
    #data.loc[:,'Sum'] = data.loc[:,['Surf1','Surf2']].sum(axis=1, min_count=1)
    dfMerge.loc[:,'Peptide_observed'] = dfMerge.loc[:,['Peptide_observed_x','Peptide_observed_y']].sum(axis=1, min_count=1)
    dfMerge.loc[:,'Total_Batches'] = dfMerge.loc[:,['Total_Batches_x','Total_Batches_y']].max(axis=1)
    
    peptideCntDict = dict(zip(dfMerge.L_peptide,dfMerge.Peptide_observed))
    batchCntDict = dict(zip(dfMerge.L_peptide,dfMerge.Total_Batches))
    return peptideCntDict, batchCntDict

def loess_newRTs(df, refDF):
    X=list(df.ReferenceRT)
    Y=list(df.NewRT)
    
    target = list(refDF.RT)

    # Build a LOESS model and calibrate RTs
    mod = rLoess(FloatVector(X), FloatVector(Y))  # LOESS model based on the shared peptides
    cal_target = rPredict(mod, FloatVector(target))  # Calibration is applied to whole peptides of the current run
    refDF["calibratedRTs"] = cal_target

    #idx = (~ref.isna()) & (~target.isna())
    np.where(refDF["calibratedRTs"].isna(),refDF["RT"], refDF["calibratedRTs"])



#This function updates the overlapped entries
#Rules: Check for the maximum Jscore of entries in Reference Library and New Library
#If Reference Library JScore >= New Library JScore
#Keep the reference entry
#Keep the new entry

#input are output of parseLib_modelDF each index separated 
def QC_newLib_Update(dfR_T, dfL_T, modelRTDF):#input are reference library dataframe (target only) and new library dataframe
    
    #OverlappedPeptide Make one dataframe 
    overlapPep = list(modelRTDF.L_peptide)
    
    overlap_R_DF = dfR_T.loc[dfR_T.L_peptide.isin(overlapPep)]
    overlap_L_DF = dfL_T.loc[dfL_T.L_peptide.isin(overlapPep)]

    
    #add this horizontally after deciding on Jscore
    keepDF_ref = overlap_R_DF[['scan', 'charge', '[M+H]+', 'precursorMZ', 'L_ID', 'L_peptide','mod_pos_line', 'RT']]
    
    #this is very important as these analyze columns are required later 
    analyzeCols = ['L_peptide','L_protein','spectrum_line', 'jscore_line', 'max_JScore', 'psm_no_line', 'm/z', 'intensity', 'Lib_Description']
    ref_check = overlap_R_DF[analyzeCols]
    new_check = overlap_L_DF[analyzeCols]

    #merge 2 dataframes based on L_peptide
    overlapMerged = ref_check.merge(new_check, on="L_peptide")
    

    #df['que'] = np.where((df['one'] >= df['two']) & (df['one'] <= df['three']), df['one'], np.nan)
    overlapMerged["Keep_Type"] = np.where(overlapMerged.max_JScore_x >= overlapMerged.max_JScore_y, "Ref","New")

    #make the column replacement dictionaries for reference and new library
    #_x is reference
    #_y is new library
    
    refColDict = {}
    newColDict = {}
    for cols in overlapMerged.columns:
        if "_x" in cols:
            colName = cols.split("_x")[0]
            refColDict[cols] = colName

        if "_y" in cols:
            colName = cols.split("_y")[0]
            newColDict[cols] = colName
            
            
    #generate dataframe to keep the columns for respective keep type
    #reference
    refKeepDF = overlapMerged.loc[overlapMerged.Keep_Type == "Ref"]
    refKeepDF2 = refKeepDF.rename(columns = refColDict)
    write_log ("Selected overlapped entries from reference library = ",refKeepDF2.shape[0])
    #new library
    newKeepDF = overlapMerged.loc[overlapMerged.Keep_Type == "New"]
    newKeepDF2 = newKeepDF.rename(columns = newColDict)
    write_log ("Selected overlapped entries from new library = ",newKeepDF2.shape[0])
    #make final matrix with only analyze columns for reference and new library
    refKeepDF2Final = refKeepDF2[analyzeCols]
    newKeepDF2Final = newKeepDF2[analyzeCols]
    
    #merge reference with new library overlapped entries based on Jscore entries
    overlapFinalDecided = refKeepDF2Final.append(newKeepDF2Final)
    
    #bring keepDF_ref and merge that with overlap matrix based on decision of Jscore 
    db_overlaped_merged = keepDF_ref.merge(overlapFinalDecided, on="L_peptide")
    
    write_log ("Percentage of overlapped peptide selected from reference library = ",np.round(refKeepDF2.shape[0]/db_overlaped_merged.shape[0]*100,2),"%")
    write_log ("Percentage of overlapped peptide selected from new library = ",np.round(newKeepDF2Final.shape[0]/db_overlaped_merged.shape[0]*100,2),"%")
    
    return db_overlaped_merged


#very simple model where we keep all reference peptide and just update the new entries RT

def QC_newLib_Update_keep_all_reference(dfR_T, dfL_T, f):#ref library; new library; and lowess curve
    #OverlappedPeptide Make one dataframe 
    
    new_peptide_Df = dfL_T[~dfL_T.L_peptide.isin(dfR_T.L_peptide)]
    new_peptide_Df["calibratedRTs"] = new_peptide_Df.apply(lambda x: np.round(f(x.RT),3), axis=1)
    
    
    return new_peptide_Df
    


'''
columns :['scan', 'charge', '[M+H]+', 'precursorMZ', 'L_ID', 'L_peptide',
       'mod_pos_line', 'spectrum_line', 'jscore_line', 'psm_no_line', 'RT',
       'm/z', 'intensity', 'Lib_Description', 'calibratedRTs']

'''

def mergeLibrary(df,specLibFolder,notes, libtypename, peptideCntDict, batchCntDict):
    mz_cols = list(df.columns)
    
    proton = 1.00727646677
    write_log ("  Generating .ms2 files\n")
    now = datetime.now()
    #print (prev_mod)
    write_log("  now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index
    
    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"
    new_ms2_file = specLibFolder+"/intermediate/jumplib_human_{}_target.splib".format(libtypename)
    # new_ms2_file = specLibFolder+"/Spectral_Library_theo_mz_exp_intensity.spLib" #modified the extension name to ms2pep
    write_log ("  Target IDs spectral library is being created", new_ms2_file)
    
    with open(new_ms2_file,"w") as new_ms2:
        new_ms2.write(header_ms2)
    
            
        #idea is to convert two list to one dictionary and keep updating dictionary as it sees another dictionary
        # keys = ['a', 'b', 'c']
        # values = [1, 2, 3]
        # dictionary = dict(zip(keys, values))
        #https://stackoverflow.com/questions/209840/how-do-i-convert-two-lists-into-a-dictionary
        #Update one dictionary from another https://stackoverflow.com/questions/10461531/merge-and-sum-of-two-dictionaries
        #x = {'both1': 1, 'both2': 2, 'only_x': 100}
        #y = {'both1': 10, 'both2': 20, 'only_y': 200}
        #print {k: x.get(k, 0) + y.get(k, 0) for k in set(x) | set(y)}
       
        
        id_no = 1
        for row in df.itertuples():
            mzIntDict = {}
            scan = row[mz_cols.index("scan")+1]
            charge = row[mz_cols.index("charge")+1]
            L_ID = row[mz_cols.index("L_ID")+1]
            L_peptide = row[mz_cols.index("L_peptide")+1]
            
            L_protein = str(row[mz_cols.index("L_protein")+1])
            
            mod_pos_line = row[mz_cols.index("mod_pos_line")+1]
            spectrum_line = row[mz_cols.index("spectrum_line")+1]
            jscore_line = row[mz_cols.index("jscore_line")+1]
            
            MH_proton = row[mz_cols.index("[M+H]+")+1]
            prec_MZ = row[mz_cols.index("precursorMZ")+1]
            
            RT = row[mz_cols.index("RT")+1]
            psm_cnt = row[mz_cols.index("psm_no_line")+1]
            mz = row[mz_cols.index("m/z")+1]
            intensity = row[mz_cols.index("intensity")+1]
            cal_RT = float(row[mz_cols.index("calibratedRTs")+1])
            libraryNotes = row[mz_cols.index("Lib_Description")+1]
            

            
            id_no_final = "p"+str(id_no).zfill(7)
            new_ms2.write("S\t"+str(id_no)+"\t"+str(id_no)+"\t"+str(prec_MZ)+"\n") #this is for fake scan for second search with comet
            new_ms2.write("Z\t"+str(charge)+"\t"+str(MH_proton)+"\n")
            if libraryNotes == notes:
                new_ms2.write("L\tID_with_Modification\t"+L_ID+"\t"+L_peptide+"\n")
            else:
                new_ms2.write("L\tID_with_Modification\t"+id_no_final+"\t"+L_peptide+"\n")
            
            new_ms2.write("L\tProtein\tRepresentative\t"+L_protein+"\n")

            for mods in mod_pos_line:
                new_ms2.write(mods+"\n")
            new_ms2.write(spectrum_line+"\n")
            new_ms2.write(jscore_line+"\n")
            new_ms2.write("L\tMH+\t"+str(MH_proton)+"\n")
            new_ms2.write("L\tPrecursor\tz="+str(charge)+"\t"+str(prec_MZ)+"\n")
            if libraryNotes != notes:
                if np.isnan(cal_RT):
                    new_ms2.write("L\tRT\t"+str(np.round(RT,3))+"\n") 
                else:
                    new_ms2.write("L\tRT\t"+str(np.round(cal_RT,3))+"\n")
            else:
                new_ms2.write("L\tRT\t"+str(np.round(RT,3))+"\n")      
            new_ms2.write(psm_cnt+"\n")
            new_ms2.write("L\tLibraryNotes"+"\t"+libraryNotes+"\n")

            #track the peptides count across batches
            totalBatchesPepCnt = str(int(peptideCntDict[L_peptide]))
            #+1 updates the batch number
            totalBatchesAnalyzed = str(int(batchCntDict[L_peptide])+1)
            new_ms2.write("L\tTotalBatches\t{}\n".format(totalBatchesAnalyzed))
            new_ms2.write("L\tPeptideLibraryCount\t{}\n".format(totalBatchesPepCnt))

            for index, val in enumerate(mz):
                new_ms2.write(str(val)+"\t"+str(intensity[index])+"\n")
            id_no+=1  ###updates ID number
            
        write_log ("  Done ...\n")



def decoySpecLibrary(df, specLibFolder, d, libtypename):
    #Decoy precursor generation  
    #d can be parameter in future

    df["Decoy_prec_mz"] = df.apply(precSwap, d=d, axis=1)

    scan_cnt = df.shape[0]
    mz_cols = list(df.columns)
    proton = 1.00727646677
    write_log ("  Generating .ms2 files\n")
    now = datetime.now()
    #print (prev_mod)
    write_log("  now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index
    
    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"
    new_ms2_file = specLibFolder+"/intermediate/jumplib_human_{}_decoy.splib".format(libtypename)
    # new_ms2_file = specLibFolder+"/Spectral_Library_Decoy.spLib" #modified the extension name to ms2pep
    write_log ("  Decoy spectral library is being created")
    
    with open(new_ms2_file,"w") as new_ms2:
        new_ms2.write(header_ms2)

        for row in df.itertuples():
            mzIntDict = {}
            scan = str(int(row[mz_cols.index("scan")+1])+scan_cnt)
            L_ID = "Decoy_"+row[mz_cols.index("L_ID")+1]
            L_peptide = "Decoy_"+row[mz_cols.index("L_peptide")+1]
            L_protein = "Decoy_"+row[mz_cols.index("L_protein")+1] 
            precursorMZ = float(row[mz_cols.index("Decoy_prec_mz")+1]) #because it is decoy precursor ions using precursor swap
            
            if "ms2_mz_int_array" in mz_cols:
                ms2_mz_int_array = row[mz_cols.index("ms2_mz_int_array")+1] #this is dictionary of mz (sorted) and intensity (summed)
            else:
                mz = row[mz_cols.index("m/z")+1]
                intensity = row[mz_cols.index("intensity")+1]
                

            charge = int(row[mz_cols.index("charge")+1])
            RT = float(row[mz_cols.index("RT")+1])

            massNeutral = (precursorMZ*charge) - ((charge-1)*proton)
            

            new_ms2.write("S\t"+scan+"\t"+scan+"\t"+str(precursorMZ)+"\n")
            new_ms2.write("Z\t"+str(charge)+"\t"+str(massNeutral)+"\n")
            new_ms2.write("L\tID_with_Modification_Decoy\t"+L_ID+"\t"+L_peptide+"\n")
            new_ms2.write("L\tProtein\tRepresentative\t"+L_protein+"\n")
            
            new_ms2.write("L\tMH+_Decoy\t"+str(massNeutral)+"\n")
            new_ms2.write("L\tPrecursor_Decoy\tz="+str(charge)+"\t"+str(precursorMZ)+"\n")
            new_ms2.write("L\tRT_Decoy\t"+str(RT)+"\n")

            if "ms2_mz_int_array" in mz_cols:
                for key,value in ms2_mz_int_array.items():
    #                 print (key,"\t",value)
                    new_ms2.write(str(key)+"\t"+str(int(value))+"\n")
            else:
                for index, val in enumerate(mz):
                    new_ms2.write(str(val)+"\t"+str(intensity[index])+"\n")

    write_log ("  Done ...\n")


# def extract_merge_L_ID_protTable(refLib, newLib, outputFolder):
#     refID_File = dirname(refLib)+"/intermediate/peptide_protein_map_library.ppml"
#     newID_File = dirname(newLib)+"/intermediate/peptide_protein_map_library.ppml"

#     df_Ref = pd.read_csv(refID_File, delimiter="\t")
#     df_New = pd.read_csv(newID_File, delimiter="\t")

#     df_merge = df_Ref.append(df_New)
#     df_merge_noDup = df_merge.drop_duplicates(subset=["Peptide","Protein Accession #"])

#     df_merge_noDup.to_csv(outputFolder+"/intermediate/peptide_protein_map_library.ppml", sep="\t", index=None)


#Implementation of precursor swap Spectra ST algorithm

# def getDecoySpectrum(targetDF, precMZ, decoy_search, charge):
#     # targetDF.loc[(targetDF.precursorMZ >= precMZ-decoy_scan) & (targetDF.precursorMZ < precMZ+decoy_scan)]
#     select_mz_DF = targetDF.loc[(targetDF.precursorMZ.between((precMZ-decoy_search),(precMZ+decoy_search))) & (targetDF.charge == charge) ]
#     check=select_mz_DF.loc[select_mz_DF.precursorMZ == np.max(select_mz_DF.precursorMZ)]
#     if check.shape[0] > 1:
#         check = check.iloc[0:1]

#     return check



