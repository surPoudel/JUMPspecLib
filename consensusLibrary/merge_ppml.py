import pandas as pd
import numpy as np
from logFunctions import *
pd.options.mode.chained_assignment = None  # default='warn'


def mergeppml(ref_ppmlfile, new_ppmlfile):
    write_log("****** Merging ppml files for the two libraries *******\n")
    #reads file to dataframe
    write_log("Reading ppml files for reference and new libraries")

    write_log("  Reading {} file for reference".format(ref_ppmlfile))
    ref_id_allDF = pd.read_csv(ref_ppmlfile, delimiter="\t")

    write_log("  Reading {} file for reference".format(new_ppmlfile))
    new_id_allDF = pd.read_csv(new_ppmlfile, delimiter="\t")
    
    write_log("  Reference library ppml all protein entries = {}".format(ref_id_allDF.shape[0]))
    write_log("  New library ppml all protein entries = {}".format(new_id_allDF.shape[0]))
    

    write_log("  Generating dictionaries")
    write_log("    1. Protein accession --> Protein Group")
    write_log("    2. Protein group num --> Protein Group")
    write_log("    3. Gene Name --> Protein group num\n")

    #extracts number from Protein_grp column 
    ref_id_allDF["Protein_grp_num"] = ref_id_allDF["Protein_grp"].str.extract('(\d+)').astype(int)
    new_id_allDF["Protein_grp_num"] = new_id_allDF["Protein_grp"].str.extract('(\d+)').astype(int)

    #define dictionary with key as protein accession and value as protein group number
    ref_protGrpDict =dict(zip(ref_id_allDF["Protein Accession #"],ref_id_allDF["Protein Group#"]))
    
    #make a dictionary of protein group number and maximum sub group value
    #The idea is to use this dictionary to solve CASE3 sub group values
    #Step 1: Drop duplicates from Protein_grp_num but keep maximum from Protein_sub_group
    grpNumSubGroupMaxRef = ref_id_allDF.groupby('Protein_grp_num')['Protein_sub_group'].max().reset_index()
    maxSubGroupRefDict = dict(zip(grpNumSubGroupMaxRef.Protein_grp_num, grpNumSubGroupMaxRef.Protein_sub_group))
    
    #similarly make GN dictionary from reference file
    #If there is a different group name in new analysis that actually have same gene name from reference, we will assign the same group name of reference for that protein
    #this is critical because there could be different group number for same protein in new database as this is data dependent
    #Steps
    #First drop NAs from GN
    ref_id_allDF_dropNA = ref_id_allDF.dropna(subset=["GN"])
    ref_GN_GroupNo =dict(zip(ref_id_allDF_dropNA["GN"],ref_id_allDF_dropNA["Protein_grp_num"]))
    
    
    #Generate groupby PeptideSeqWithRealDelMass and protein group for reference and new library
    #x = ref y =new
    mergeRefNew=ref_id_allDF.merge(new_id_allDF, how="outer", on=["PeptideSeqWithRealDelMass","Protein Accession #"])
    
    write_log("  A total of {} ppml entries are under evaluation".format(mergeRefNew.shape[0]))
    #returns merged dataframe and 3 dictionaries
    return mergeRefNew, ref_protGrpDict, maxSubGroupRefDict, ref_GN_GroupNo



def updateSubGroup(new_id_case3NRFinal,maxSubGroupRefDict):

    #use the dictionary maxSubGroupRefDict to update the sub group value

    grp_num_list = list(new_id_case3NRFinal.Protein_grp_num_x)
    new_sub_group_list = []
    update = 0
    for x in range(0,len(grp_num_list)-1):
        grpVal = grp_num_list[x]
        if grpVal == grp_num_list[x+1]:
            update +=1
            refMaxValue = maxSubGroupRefDict[grpVal]
            new_sub_group_list.append(refMaxValue+update)
        else:
            update = 1
            refMaxValue = maxSubGroupRefDict[grpVal]
            new_sub_group_list.append(refMaxValue+update)
    grpVal = grp_num_list[-1]
    refMaxValue = maxSubGroupRefDict[grpVal]
    new_sub_group_list.append(refMaxValue+update)
    
    return new_sub_group_list



def updateProteinGroupNo(new_id_case4,mergeRefNew):
    #For remaining group number in New search, update group number only (maximum reference group number + 1) and keep the same sub group number

    #Steps 1. first make a dataframe that only has proteins without group names in reference from newProtNewPepDF_sorted
    #Step 2. Create the replace group number value. The group number now is updated based on maximum reference group
    #Step 3. Keep the same subtgroup as the new database results

    ########making a dictionary
    replaceGroupNoDict = {}
    sortedNewGroupNoList = list(new_id_case4.Protein_grp_num_y)
    maxProtGrp = np.max(mergeRefNew.Protein_grp_num_x)
    n=1

    for key in sortedNewGroupNoList:
        if key not in replaceGroupNoDict.keys():
            newvalue = maxProtGrp+n
            replaceGroupNoDict[key]=newvalue
            n+=1
    ############
    
    return replaceGroupNoDict



def ppmlPeptideEvaluation(ref_ppmlfile, new_ppmlfile):
    #mergeppml function is called 
    write_log("Two ppml dataframe are horizontally merged (outer merge). Evaluating different cases of new peptide to get the combined ppml file")

    mergeRefNew, ref_protGrpDict, maxSubGroupRefDict, ref_GN_GroupNo = mergeppml(ref_ppmlfile, new_ppmlfile)
    
    #New peptide 
    #Belong to old protein (old DB group and old DB subgroup)
    #Belong to new protein (new DB group number and new subgroup) – Rename new DB group number
    #Belongs to both (old DB group number and add new DB subgroup no)
    
    #extract dataframe that have new peptides
    write_log("  Keeping Reference ppml file that have Protein Group#")
    old_id_case1_keep = mergeRefNew.loc[~mergeRefNew["Protein Group#_x"].isnull()] #we don't need to do anything with this group
    write_log("  Case 1 evaluated. Total of {} ppml entries are retained.".format(old_id_case1_keep.shape[0]))
    new_id_case2_map = mergeRefNew.loc[mergeRefNew["Protein Group#_x"].isnull()] #we need to further analyze this group
    write_log("Total of {} ppml entries are under further evaluation.".format(new_id_case2_map.shape[0]))


    #New peptide that belongs to old protein
    #no need to worry about this as the protein group and subgroups are assigned
    #CASE2 analysis
    #create a new workign df
    # new_id_case2_map = new_id_analyze.copy(deep=False)
    new_id_case2_map["Protein Group#_x"] = new_id_case2_map["Protein Accession #"].map(ref_protGrpDict)
    write_log("  The ppml new peptides entries that have Protein Accession # are now mapped to Protein Group#")
    #separate the matched population we do not need to do anything
    new_id_case2 = new_id_case2_map.loc[~new_id_case2_map["Protein Group#_x"].isnull()]
    write_log("  Case 2 evaluated. Total of {} ppml entries are retained.".format(new_id_case2.shape[0]))
    
    new_id_case3_analysis = new_id_case2_map.loc[new_id_case2_map["Protein Group#_x"].isnull()]
    write_log("Total of {} ppml entries are under further evaluation.".format(new_id_case3_analysis.shape[0]))


    #Step 1: sort the dataframe first by Protein_grp_num_y and Protein_sub_group_y
    write_log("  Sort the dataframe first by Protein_grp_num and Protein_sub_group")
    new_id_case3_sorted = new_id_case3_analysis.sort_values(['Protein_grp_num_y', 'Protein_sub_group_y'], ascending=[True,True])
    write_log("  The ppml new peptides entries that have Gene Name are now mapped to Protein Group number")
    
    #Now map the GN of new data with the group no from reference data
    new_id_case3_sorted["Protein_grp_num_x"] = new_id_case3_sorted.GN_y.map(ref_GN_GroupNo)
    #CASE 3
    new_id_case3 = new_id_case3_sorted.loc[~new_id_case3_sorted.Protein_grp_num_x.isnull()]
    
    #drop duplicates by protein accession
    new_id_case3NRFinal = new_id_case3.drop_duplicates(subset=["Protein Accession #"], keep="first")
    write_log("  The ppml entries with unique Protein Accession are kept. Maximum sub group value are assigned using maxSubGroupRefDict")

    
    #update sub group and returns the list: Input = new_id_case3NRFinal
    new_sub_group_list = updateSubGroup(new_id_case3NRFinal,maxSubGroupRefDict)
    write_log("  Sub group value are updated for each Protein_grp following maximum sub group number")
   
    # new_id_case3NRFinal = new_id_case3NR.copy(deep=False)
    new_id_case3NRFinal["update_sub_group"] = new_sub_group_list
    #Update Protein Group#_x
    #'{0:03}'.format(1)

    new_id_case3NRFinal["update_sub_group"] = new_id_case3NRFinal["update_sub_group"].apply(lambda x: "{0:03}".format(x))
    new_id_case3NRFinal["Protein_grp_num_x"] = new_id_case3NRFinal["Protein_grp_num_x"].apply(lambda x: str(int(x)).zfill(7))
    new_id_case3NRFinal["Protein Group#_x"] = "SJPG"+new_id_case3NRFinal["Protein_grp_num_x"]+"."+new_id_case3NRFinal["update_sub_group"].astype("str")
    write_log("  Protein Group# is generated using Protein Group num (GN) and updated sub group")
    #make a dictionary of "Protein Accession #" and "Protein Group#_x" and use it back to map new_id_case3 DF
    gn_derived_protGrp_accesion_dict = dict(zip(new_id_case3NRFinal["Protein Accession #"],new_id_case3NRFinal["Protein Group#_x"]))
    write_log("  A dictionary that has Protein Accession# --> Protein Group# is generated")
    
    # new_id_case3Final = new_id_case3.copy(deep=False)
    new_id_case3["Protein Group#_x"] = new_id_case3["Protein Accession #"].map(gn_derived_protGrp_accesion_dict)
    write_log("  Protein Accession# from new ppml entries are updated with Protein Group# using this dictionary ")

    write_log("  Case 3 ppml entried that have Gene Name are evaluated. Total of {} ppml entries are retained.".format(new_id_case3.shape[0]))

    #all null reference new_id_case3_sorted
    #these are the peptides that are only present in new database mapped to new proteins
    #idea is to give New SJPG number to this protein
    #The reference number will always be constant, we will only update new additions
    #CASE 2
    #Belong to new protein (new DB group number and new subgroup) – Rename new DB group number

    #CASE 4
    new_id_case4 = new_id_case3_sorted.loc[new_id_case3_sorted.Protein_grp_num_x.isnull()]
    write_log("Total of {} ppml entries are under further evaluation.".format(new_id_case4.shape[0]))

    #update new group number after the maximum reference Group number
    replaceGroupNoDict = updateProteinGroupNo(new_id_case4, mergeRefNew)
    
    write_log("  Protein group number are updated for these ppml entries. The reference group number for these entries starts after the maximum Protein Group number")
    
    # new_id_case4_solved = new_id_case4.copy(deep=False)
    new_id_case4["Protein_grp_num_x"] = new_id_case4.Protein_grp_num_y.map(replaceGroupNoDict)
    new_id_case4["Protein_sub_group_x"] = new_id_case4.Protein_sub_group_y

    new_id_case4["Protein_sub_group_x"] = new_id_case4["Protein_sub_group_x"].apply(lambda x: "{0:03}".format(int(x)))
    new_id_case4["Protein_grp_num_x"] = new_id_case4["Protein_grp_num_x"].apply(lambda x: "{0:07}".format(int(x)))
    new_id_case4["Protein Group#_x"] = "SJPG"+new_id_case4["Protein_grp_num_x"]+"."+new_id_case4["Protein_sub_group_x"].astype("str")
    
    write_log("  Sub group level of protein are borrowed from new ppml entries. The referencr Protein Group# is updated based on Protein group num and sub group number")

    write_log("All levels of ppml entries evaluated and new protein groups are updated based on Reference and New ppml\n")
    
    return old_id_case1_keep, new_id_case2, new_id_case3, new_id_case4



def reformatDF(df, colslist, replaceDict):
    #Main DF cols Important notes on column structures of all cases dataframes
    '''
    Index(['Peptides', 'PeptideSeqWithRealDelMass', 'Protein Group#',
           'Protein Accession #', 'Protein Description', 'GN', 'Fate',
           'Protein_grp', 'Protein_sub_group'],
          dtype='object')

    look at all the DFs and clean them
    have peptide information 
    keep cols
    ['Peptides_x', 'PeptideSeqWithRealDelMass', 'Protein Group#_x',
          'Protein Accession #', 'Protein Description_x', 'GN_x', 'Fate_x',
          'Protein_grp_x', 'Protein_sub_group_x']
    old_id_case1_keep.columns

    no peptide in reference 
    so use new library peptide information
    use protein group to make 'Protein_grp_x','Protein_sub_group_x'
    ['Peptides_y','PeptideSeqWithRealDelMass' ,'Protein Group#_x', 'Protein Accession #','Protein Description_y', 'GN_y',
          'Fate_y', 'Protein_grp_x','Protein_sub_group_x']
    new_id_case2.columns

    use protein group to make 'Protein_grp_x','Protein_sub_group_x'
    ['Peptides_y','PeptideSeqWithRealDelMass' ,'Protein Group#_x', 'Protein Accession #','Protein Description_y', 'GN_y',
          'Fate_y', 'Protein_grp_x','Protein_sub_group_x']
    new_id_case3.columns

    #use protein group to make 'Protein_grp_x','Protein_sub_group_x'
    ['Peptides_y','PeptideSeqWithRealDelMass' ,'Protein Group#_x', 'Protein Accession #','Protein Description_y', 'GN_y',
          'Fate_y', 'Protein_grp_x','Protein_sub_group_x']
    new_id_case4.columns

    '''


    df_Final = df[colslist]
    df_Final.rename(columns=replaceDict, inplace=True)
    return df_Final



def consensus_ppml(ref_ppmlfile, new_ppmlfile):
    
    old_id_case1_keep,new_id_case2,new_id_case3,new_id_case4=ppmlPeptideEvaluation(ref_ppmlfile, new_ppmlfile)
    
    '''
    Good for 
    new_id_case2
    new_id_case3Final
    new_id_case4_solved
    '''
    requiredMainCols = ['Peptides', 'Protein Group#',
            'Protein Description', 'GN', 'Fate']
    goodcols_2use= ['Peptides_y','Protein Group#_x', 
                    'Protein Description_y', 'GN_y','Fate_y']


    colReplaceDict = {}
    for i,cols in enumerate(goodcols_2use):
        colReplaceDict[cols]=requiredMainCols[i]



    #For old IDs
    old_peptide_cols = {'Peptides_x':'Peptides','Protein Group#_x':'Protein Group#',
                       'Protein Description_x':'Protein Description','GN_x':'GN','Fate_x':'Fate'}


    reqdCols = ['Peptides_x', 'PeptideSeqWithRealDelMass', 'Protein Group#_x',
          'Protein Accession #', 'Protein Description_x', 'GN_x', 'Fate_x']


    '''For:
    new_id_case2
    new_id_case3
    new_id_case4
    ''' 
    reqdColsNew = ['Peptides_y','PeptideSeqWithRealDelMass' ,'Protein Group#_x', 
                   'Protein Accession #','Protein Description_y', 'GN_y','Fate_y']


    old_id_case1_keepFinal = reformatDF(old_id_case1_keep, reqdCols, old_peptide_cols)
    new_id_case2Final = reformatDF(new_id_case2, reqdColsNew, colReplaceDict)
    new_id_case3Final = reformatDF(new_id_case3, reqdColsNew, colReplaceDict)
    new_id_case4Final = reformatDF(new_id_case4, reqdColsNew, colReplaceDict)

    final_ppml_DF =old_id_case1_keepFinal.append(new_id_case2Final.append(new_id_case3Final.append(new_id_case4Final)))
    
    write_log("All cases specific dataframes are now merged together following the proper reformating for final ppml display.")

    return final_ppml_DF



def gen_merged_ppml(ref_ppmlfile, new_ppmlfile, specLibFolder):
    
    df = consensus_ppml(ref_ppmlfile, new_ppmlfile)

    
    df[["Protein_grp","Protein_sub_group"]] = df["Protein Group#"].str.split(".",expand=True)
    df["Protein_sub_group"]=df.Protein_sub_group.astype("int")
    dfNR = df.loc[df.groupby(['PeptideSeqWithRealDelMass','Protein_grp']).Protein_sub_group.idxmin()]
    #Generate groupby PeptideSeqWithRealDelMass and protein group 
    dfNR["Protein_grp_num"] = dfNR["Protein_grp"].str.extract('(\d+)').astype(int)
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
    
    makedirectory(specLibFolder+"/intermediate")
    
    df.to_csv(specLibFolder+"/intermediate/id_all_pep.ppml", sep="\t", index=None)
    dfNR2.to_csv(specLibFolder+"/intermediate/id_uni_pep.ppml",sep="\t", index=None)

    write_log("The unique protein level and all protein level ppml files are generated and stored as {} {}".format(specLibFolder+"/intermediate/id_all_pep.ppml",specLibFolder+"/intermediate/id_uni_pep.ppml"))


