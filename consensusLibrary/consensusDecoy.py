import pandas as pd
import numpy as np
import os, sys, glob, re
from datetime import datetime
from logFunctions import *


def precSwap(row, d): #d is the distance in dalton #df is the target dataframe 
    scan = int(row.scan)
    precmz = float(row.precursorMZ)
    charge = int(row.charge)
    if scan%2 == 0:
        precmz_decoy = float(row.precursorMZ)+d
    else:
        precmz_decoy = float(row.precursorMZ)-d
    # precmz_Min = float(row.precursorMZ)-d
    
    #subset the dataframe to get the precMZ within the range of max and min (that is maximum), also need to check that it is not same precMZ as target
    #df2 = df.loc[(df.precursorMZ.between(precmz_Min, precmz_Max)) & (df.charge == charge) & (df.precursorMZ != precmz)]
    #decoyMZ =  np.max(df2.precursorMZ) #creating the precursor mass with maximum m/z within the range
    return precmz_decoy


def decoySpecLibrary(df_all, specLibFolder, d, libtypename):
    #Decoy precursor generation  
    #d can be parameter in future
    #some of the psms in ID.txt are not present in id_all_pep.txt which is used to make ppml files
    #to avoid such problems, we remove all the entries that do not have Protein columns

    df = df_all.dropna(subset=["Protein"])

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
    new_ms2_file = specLibFolder+"/intermediate/jumplib_human_{}_decoy.splib".format(libtypename) #modified the extension name to ms2pep
    write_log ("  Decoy spectral library is being created")
    
    with open(new_ms2_file,"w") as new_ms2:
        new_ms2.write(header_ms2)

        for row in df.itertuples():
            mzIntDict = {}
            scan = str(int(row[mz_cols.index("scan")+1])+scan_cnt)
            L_ID = row[mz_cols.index("L_ID")+1]
            peptide_ID = row[mz_cols.index("Peptide_ID")+1]
            precursorMZ = float(row[mz_cols.index("Decoy_prec_mz")+1]) #because it is decoy precursor ions using precursor swap
            ms2_mz_int_array = row[mz_cols.index("ms2_mz_int_array")+1] #this is dictionary of mz (sorted) and intensity (summed)
            charge = int(row[mz_cols.index("charge")+1])
            RT = float(row[mz_cols.index("RT")+1])
            protein = row[mz_cols.index("Protein")+1]
            massNeutral = (precursorMZ*charge) - ((charge-1)*proton)
            

            new_ms2.write("S\t"+scan+"\t"+scan+"\t"+str(precursorMZ)+"\n")
            new_ms2.write("Z\t"+str(charge)+"\t"+str(massNeutral)+"\n")
            new_ms2.write("L\ttID_with_Modification_Decoy\tDecoy_"+L_ID+"\tDecoy_"+peptide_ID+"\n")
            new_ms2.write("L\tProtein\tRepresentative\tDecoy_"+protein+"\n")
            new_ms2.write("L\tMH+_Decoy\t"+str(massNeutral)+"\n")
            new_ms2.write("L\tPrecursor_Decoy\tz="+str(charge)+"\t"+str(precursorMZ)+"\n")
            new_ms2.write("L\tRT_Decoy\t"+str(RT)+"\n")

            for key,value in ms2_mz_int_array.items():
#                 print (key,"\t",value)
                new_ms2.write(str(key)+"\t"+str(int(value))+"\n")
    write_log ("  Done ...\n")





def getDecoySpectrum(targetDF, precMZ, decoy_search, charge):
    # targetDF.loc[(targetDF.precursorMZ >= precMZ-decoy_scan) & (targetDF.precursorMZ < precMZ+decoy_scan)]
    select_mz_DF = targetDF.loc[(targetDF.precursorMZ.between((precMZ-decoy_search),(precMZ+decoy_search))) & (targetDF.charge == charge) ]

    #check if the select_mz_DF has only one row. If so, it is the precursor scan
    totalCandidates = select_mz_DF.shape[0]
    updateDecoy = 1
    while totalCandidates == 1:
        updateDecoy +=1
        decoy_search = decoy_search*updateDecoy
        select_mz_DF = targetDF.loc[(targetDF.precursorMZ.between((precMZ-decoy_search),(precMZ+decoy_search))) & (targetDF.charge == charge) ]
        totalCandidates = select_mz_DF.shape[0]

    select_mz_DF["abs_prec_diff"] = (select_mz_DF.precursorMZ - decoy_search).abs()
    check=select_mz_DF.loc[select_mz_DF.precursorMZ == np.max(select_mz_DF.abs_prec_diff)]
    if check.shape[0] > 1:
        check = check.iloc[0:1]

    return check

def getDecoySpectrum_SpectraST(targetDF, precMZ, decoy_search, charge, exclusion_list): #dyanmic swap_d is the range 8 Da + 1 or 2 or 3 where we look for swap precursor
    
    #This function uses precMZ to scan over all library outside the +/-decoy search range but closest precMZ value
    #the precursor MZ is not between the range of decoy_search; Δ m/z > decoy_search. We want to have same charge. Now this will give us 100s of candidates
    
    # reqd_cols = ["precursorMZ","charge","ms2_mz_int_array","intensity","scan","L_ID","Peptide_ID","Protein","RT"]
    #select by charge
    select_mz_DF_11 = targetDF.loc[targetDF["charge"] == charge]
    select_mz_DF_11["scan"]=select_mz_DF_11.scan.astype("int")
    select_mz_DF_11.set_index("scan", inplace=True)
    select_mz_DF_1 = select_mz_DF_11.loc[~select_mz_DF_11.index.isin(exclusion_list)]
    select_mz_DF_1["abs_delMZ"] = abs(np.array(select_mz_DF_1.precursorMZ)-precMZ)
    #apply cutoff |8Da|
    #Adding 200 Dalton as cutoff should help us to find at least 1 decoy candidate
    

    
    if select_mz_DF_1.shape[0] <2:
        return None
    
    else:
    
        select_mz_DF=select_mz_DF_1.loc[select_mz_DF_1["abs_delMZ"].between((decoy_search),(decoy_search+10000))] #+/- 1 dalton

#     #     #make sure we have at least one decoy found
#         attempts = 1
#         while (select_mz_DF.shape[0] < 1) and (attempts < 3):
#             # print ("Decoy not available with d = {} for precursor mz {} and charge {}".format(decoy_search,precMZ,charge))
#             decoy_search +=1
#             attempts+=1
#             select_mz_DF=select_mz_DF_1.loc[select_mz_DF_1["abs_delMZ"].between((decoy_search),(decoy_search+dynamic_swap_d))]
#             # print ("Searching precursor mz {} decoy with attempt {} and decoy mass d= {} Dalton ".format(precMZ,attempts,decoy_search))

        if select_mz_DF.shape[0] == 0:
            return None
        else:

            #sort by ascending order of absolute difference 
            select_mz_DF2 = select_mz_DF.copy().sort_values(by=["abs_delMZ"], ascending=True)
            
                        
            #select the first row
            check = select_mz_DF2.iloc[0:1]
            check.reset_index(inplace=True)
            mz_cols = list(check.columns)
            arr = check.to_numpy()
            row = arr[0]
            
            if "ms2_mz_int_array" in mz_cols:
                mz = row[mz_cols.index("ms2_mz_int_array")]
            else:
                mz_array = row[mz_cols.index("m/z")]
                intensity = row[mz_cols.index("intensity")]
                mz = dict(zip(mz_array,intensity))
            
#             mz = row[mz_cols.index("ms2_mz_int_array")]
            # intensity = list(row[mz_cols.index("intensity")])
            scan_pair = int(row[mz_cols.index("scan")])
            precMZ_pair = float(row[mz_cols.index("precursorMZ")])
#             L_ID_pair = "Decoy_"+row[mz_cols.index("L_ID")]
#             L_peptide_pair = "Decoy_"+row[mz_cols.index("Peptide_ID")]
#             L_protein_pair = "Decoy_"+row[mz_cols.index("Protein")]
            
            L_ID_pair = "Decoy_"+row[mz_cols.index("L_ID")]
            if "Peptide_ID" in mz_cols:
                L_peptide_pair = "Decoy_"+row[mz_cols.index("Peptide_ID")]
            else:
                L_peptide_pair = "Decoy_"+row[mz_cols.index("L_peptide")]
            if "Protein" in mz_cols: 
                L_protein_pair = "Decoy_"+row[mz_cols.index("Protein")] 
            else:
                L_protein_pair = "Decoy_"+row[mz_cols.index("L_protein")]
            
            RT_pair = float(row[mz_cols.index("RT")])

            # print ("for {} precursor the pair scan is {} of prec m/z {}".format(precMZ,scan_pair,precMZ_pair))

            attr_list = [mz,scan_pair,precMZ_pair,L_ID_pair,L_peptide_pair,L_protein_pair,RT_pair]
            
            return attr_list


def decoySpecLibrary_Prec_Swap_New(df, specLibFolder, d, libtypename, scan_cnt, exclusion_list, rescue_scans_list):
    #rescue_scans_list collects the list of scans (entries) that did not find the pair from that bins. These entries are 
    #collected from all bins together and they search on the whole database for their counter decoy -- this could
    #be some outlier for our decoy generation as already formed pair could be involved in this pairing 
    mz_cols = list(df.columns)
    tracker = 0
    np_arr = df.to_numpy()
    proton = 1.00727646677
    decoy_dict = {}
    for row in np_arr:
        scan = int(row[mz_cols.index("scan")])

        decoy_scan = scan+scan_cnt

        if scan not in exclusion_list:

            charge = int(row[mz_cols.index("charge")])
            precMZ = float(row[mz_cols.index("precursorMZ")])
#             L_ID = "Decoy_"+row[mz_cols.index("L_ID")]
#             L_peptide = "Decoy_"+row[mz_cols.index("Peptide_ID")]
#             L_protein = "Decoy_"+row[mz_cols.index("Protein")] 
            
            L_ID = "Decoy_"+row[mz_cols.index("L_ID")]
            if "Peptide_ID" in mz_cols:
                L_peptide = "Decoy_"+row[mz_cols.index("Peptide_ID")]
            else:
                L_peptide = "Decoy_"+row[mz_cols.index("L_peptide")]
            if "Protein" in mz_cols: 
                L_protein = "Decoy_"+row[mz_cols.index("Protein")] 
            else:
                L_protein = "Decoy_"+row[mz_cols.index("L_protein")] 

            attrlist = getDecoySpectrum_SpectraST(df, precMZ, d, charge, exclusion_list)
            if attrlist == None:
                rescue_scans_list.append(scan)
                pass
            else:
                
                
                mz,scan_pair,precMZ_pair,L_ID_pair,L_peptide_pair,L_protein_pair,RT_pair = attrlist
                decoy_scan_pair = str(scan_pair+scan_cnt) 
                
                key1 = "{}_{}".format(scan_pair,precMZ_pair)

                exclusion_list.append(scan)
                exclusion_list.append(scan_pair)

                RT = float(row[mz_cols.index("RT")])
                massNeutral = (precMZ*charge) - ((charge-1)*proton)
                LibraryNotes = "Target Pair = {}; prec mz = {}".format(scan_pair, precMZ_pair)          
                LibraryNotes_pair = "Target Pair = {}; prec mz = {}".format(scan, precMZ)          


#                 mz_pair = row[mz_cols.index("ms2_mz_int_array")]
                
                if "ms2_mz_int_array" in mz_cols:
                    mz_pair = row[mz_cols.index("ms2_mz_int_array")]
                else:
                    mz_array = row[mz_cols.index("m/z")]
                    intensity = row[mz_cols.index("intensity")]
                    mz_pair = dict(zip(mz_array,intensity))


                massNeutral_pair = (precMZ_pair*charge) - ((charge-1)*proton)
                
                decoy_dict[key1] = [decoy_scan, precMZ, charge, massNeutral, L_ID, L_peptide,L_protein, RT,mz,LibraryNotes]
                
                key2 = "{}_{}".format(scan,precMZ)
                decoy_dict[key2] = [decoy_scan_pair, precMZ_pair, charge, massNeutral_pair, L_ID_pair, L_peptide_pair,L_protein_pair, RT_pair,mz_pair,LibraryNotes_pair]
                tracker+=1

                if np.mod((tracker),10000) == 0:
                    print ("Total scans completed {}".format(tracker))

    result = [exclusion_list, rescue_scans_list, decoy_dict]     
    return result


def rescue_scan_decoy(df, rescue_scans_list, distanceDecoy):
    #rescue_scans_list collects the list of scans (entries) that did not find the pair from that bins. These entries are 
    #collected from all bins together and they search on the whole database for their counter decoy -- this could
    #be some outlier for our decoy generation as already formed pair could be involved in this pairing 

    #make a new dataframe with the scans to be rescued
    df_Rescue = df.set_index("scan")
    df_Rescue_final = df_Rescue.loc[df_Rescue.index.astype("int").isin(rescue_scans_list)]
    df_Rescue_final.reset_index(inplace=True)
    mz_cols = list(df_Rescue_final.columns)
    np_arr = df_Rescue_final.to_numpy()
    decoy_add_scan = df.shape[0]
    proton = 1.00727646677  
    decoy_dict = {}
    for row in np_arr:
        if "ms2_mz_int_array" in mz_cols:
            mz_pair = row[mz_cols.index("ms2_mz_int_array")]
        else:
            mz_array = row[mz_cols.index("m/z")]
            intensity = row[mz_cols.index("intensity")]
            mz_pair = dict(zip(mz_array,intensity))
            
#         mz_pair = row[mz_cols.index("ms2_mz_int_array")]
        # intensity = list(row[mz_cols.index("intensity")])
        scan = int(row[mz_cols.index("scan")])
        
        scan_pair = scan+decoy_add_scan
        charge = int(row[mz_cols.index("charge")])
        
        precMZ = float(row[mz_cols.index("precursorMZ")])
        precMZ_pair = precMZ+distanceDecoy #8Dalton addition
        
#         L_ID_pair = "Decoy_"+row[mz_cols.index("L_ID")]
#         L_peptide_pair = "Decoy_"+row[mz_cols.index("Peptide_ID")]
#         L_protein_pair = "Decoy_"+row[mz_cols.index("Protein")]
        
        L_ID_pair = "Decoy_"+row[mz_cols.index("L_ID")]
 

        if "Peptide_ID" in mz_cols:
            L_peptide_pair = "Decoy_"+row[mz_cols.index("Peptide_ID")]
        else:
            L_peptide_pair = "Decoy_"+row[mz_cols.index("L_peptide")]
        if "Protein" in mz_cols: 
            L_protein_pair = "Decoy_"+row[mz_cols.index("Protein")] 
        else:
            L_protein_pair = "Decoy_"+row[mz_cols.index("L_protein")] 
        
        RT_pair = float(row[mz_cols.index("RT")])

         #new_ms2.write("L\tLibraryNotes\tTarget Pair = {}; prec mz = {}\n".format(scan, precMZ))
        LibraryNotes = "Target Pair = {}; prec mz = {}".format(scan, precMZ)          
        # print ("for {} precursor the pair scan is {} of prec m/z {}".format(precMZ,scan_pair,precMZ_pair))

        massNeutral_pair = (precMZ_pair*charge) - ((charge-1)*proton) 
        
        key2 = "{}_{}".format(scan,precMZ)
        decoy_dict[key2] = [scan_pair, precMZ_pair, charge, massNeutral_pair, L_ID_pair, L_peptide_pair,L_protein_pair, RT_pair,mz_pair,LibraryNotes]
               
    return decoy_dict
       
       

def write_decoy_library(decoy_master_list, specLibFolder, d, libtypename):
    write_log ("  Generating .ms2 files\n")
    now = datetime.now()
    #print (prev_mod)
    write_log("  now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index
    
    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"
#     new_ms2_file = specLibFolder+"/intermediate/jumplib_human_{}_decoy.splib".format(cnt)
    new_ms2_file = specLibFolder+"/intermediate/jumplib_human_{}_decoy.splib".format(libtypename)
    print (new_ms2_file)
    write_log ("  Decoy spectral library is being created")
    
    tracker = 0
    
    with open(new_ms2_file,"w") as new_ms2:
        new_ms2.write(header_ms2)
        for decoy_dict in decoy_master_list:
            for key, value in decoy_dict.items():
                tracker+=1
                #decoy_master_dict[key2] = [scan_pair, precMZ_pair, charge, massNeutral_pair, L_ID_pair, L_peptide_pair,L_protein_pair, RT_pair,mz_pair]
                decoy_scan_pair, precMZ_pair, charge, massNeutral_pair, L_ID_pair, L_peptide_pair,L_protein_pair, RT_pair,mz_pair,LibraryNotes = value

                new_ms2.write("S\t{}\t{}\t{}\n".format(decoy_scan_pair, decoy_scan_pair, precMZ_pair))
                new_ms2.write("Z\t{}\t{}\n".format(charge, massNeutral_pair))

                new_ms2.write("L\tID_with_Modification_Decoy\t"+L_ID_pair+"\t"+L_peptide_pair+"\n")
                new_ms2.write("L\tProtein_Decoy\tRepresentative\t"+L_protein_pair+"\n")

                new_ms2.write("L\tMH+_Decoy\t"+str(massNeutral_pair)+"\n")
                new_ms2.write("L\tPrecursor_Decoy\tz="+str(charge)+"\t"+str(precMZ_pair)+"\n")
                new_ms2.write("L\tRT_Decoy\t"+str(RT_pair)+"\n")
                new_ms2.write("L\tLibraryNotes\t{}\n".format(LibraryNotes))

                
                for key,value in mz_pair.items():
                    new_ms2.write(str(key)+"\t"+str(int(value))+"\n")


    #                     print (scan,scan_pair,exclusion_list)
                if np.mod((tracker),10000) == 0:
                    print ("Total scans completed {}".format(tracker))

                    

    write_log ("  Done ...\n")
