import pandas as pd
import numpy as np
import pyteomics
from collections import Counter
from pyteomics import mass
from logFunctions import *
from TargetDecoy import *

def my_agg(x):
    names = {
        'Total Peptide#': x['Scan#'].count(),
        'Unique Peptide#': x['unique'].sum(),
        'PSM#': x['PSM#'].sum()}
    return pd.Series(names)


def binning(df, psm_perBin=1000,binCol = "JDscore"):
    df2 = df.copy()
    bins = int(df2.shape[0]/psm_perBin)
    df2['mzBin'] = pd.qcut(df2[binCol], q=bins)
    return df2

def bin_fixed(df, bins, bincol = "JDscore"): #bins is the list of bins
    df['bins_{}'.format(bincol)] = pd.cut(x=df[bincol], bins=bins)
    return df

def calcFDR_df(df):
    decoy = df.loc[df.Peptide.str.contains("Decoy")]
    if df.shape[0] == decoy.shape[0]:
        fdr = 100000
    else:
        fdr = decoy.shape[0]/(df.shape[0]-decoy.shape[0])*100
    write_log("  Target = {}, Decoy = {} ".format((df.shape[0]-decoy.shape[0]),decoy.shape[0]))
    return fdr


def jumpl_RT_filter(binDF,FDR):
    
    if initFDR < FDR:
        return binDF,initFDR,"No Iteration"
    else:
        try:
            maxDelRTDecoy = int(np.max(binDF.loc[binDF.Peptide.str.contains("Decoy")].deltaRT.abs()))
        except:
            maxDelRTDecoy =1000

        df_fdr_bins = binDF.loc[binDF.deltaRT.abs() < maxDelRTDecoy]
        fdr = calcFDR_df(df_fdr_bins)
        return df_fdr_bins,fdr,maxDelRTDecoy


def countTargetDecoy(df):
    target = df.loc[~df.Peptide.str.contains("Decoy")]
    decoy = df.loc[df.Peptide.str.contains("Decoy")]
    targetCnt = target.shape[0]
    decoyCnt = decoy.shape[0]
    return [targetCnt,decoyCnt]


###### OLDER function #### SD based RT tolerance
def jump_filter_library(df):
    SD=1
    #first extract data into 2 parts
    
    #highest score > 0.95 
    df2= df.copy()
    dfHighScore = df2.loc[df2.JDscore > 0.95]
    
    dfLowScore = df2.loc[df2.JDscore <= 0.95]
    #use RT filter here
    
    rt_tolerance = SD*np.std(df2.deltaRT)
    dfLowScore_rtFiltered = dfLowScore.loc[dfLowScore.deltaRT.abs() < rt_tolerance]
    
    #concatenate highest score df and rt filtered lower score df
    
    df_RT_filt = dfHighScore.append(dfLowScore_rtFiltered)
    return df_RT_filt    



def peptideLevelFiltering(psms_DF):
    #first make a column run that all the path
    psms_DF[["Run","Scan#","ppi_rank","z"]]=psms_DF["Outfile"].str.split(".",expand=True)
    #just add fraction name in Run#
    psms_DF["Run#"] = psms_DF["Run"].apply(lambda x: x.split("/")[-1])
    #drop Run column
    psms_DF.drop(columns=["Run"], inplace=True)
    
    psms_DF["m/z"] = psms_DF["measuredMH"]
    
    #unique psms level
    psms_DF_unique = psms_DF[psms_DF.Protein == psms_DF.Unique_Protein]
    
    #picks highest JDscore peptide to make a unique peptide matrix
    #For this L_ID is used. We cannot use Peptide_ID as Peptide_ID of decoy is Decoy only so we cannot group them
    #We also need to add the count of PSM# for each peptide

    write_log("  The psms level filtering is complete. The unique psms level publication file is generated. Preparing peptide level report.")
    df_peptide = psms_DF_unique.loc[psms_DF_unique.groupby('Peptide').JDscore.idxmax()].sort_values(by=["JDscore"], ascending=False)
    df_peptide_all = psms_DF.loc[psms_DF.groupby(['Peptide','Protein']).JDscore.idxmax()].sort_values(by=["JDscore"], ascending=False)
    
    
    #add PSM# count for each L_ID
    psms_count=psms_DF_unique.groupby('Peptide')['Outfile'].count().reset_index()
    psms_count2 =psms_count.rename(columns={"Outfile":"PSM#"})
    df_peptide_count = df_peptide.merge(psms_count2,how="outer",on="Peptide")
    #compute peptide FDR
    peptideFDR = calcFDR_df(df_peptide_count)
    
    df_peptide_count_all = df_peptide_all.merge(psms_count2,how="outer",on="Peptide")
    
    df_peptide_count["Peptides"] = df_peptide_count["Peptide"]
    df_peptide_count_all["Peptides"] = df_peptide_count_all["Peptide"]
    
    
    
    write_log("  The final peptide FDR = {}\n".format(peptideFDR))
    return df_peptide_count, df_peptide_count_all


def proteinLevelFiltering(peptideDF, uniqProtAllProtDict):
    
    #picks highest JDscore peptide to make a unique protein matrix
    #For this ProteinAccession is used. We cannot use Peptide_ID as Peptide_ID of decoy is Decoy only so we cannot group them
    #so we leveraged L_ID decoy information and replace the NaN in Protein Accession in the above steps
    #We also need to add the count of PSM# for each peptide
    write_log("\nThe unique peptide level report is ready. Preparing unique protein level report.")
    
    #replace unique peptide columnn with NaN for 0 so that during count it will be correct
    peptideDF["unique"]=peptideDF.unique.astype("float")
    
    df_protein = peptideDF.loc[peptideDF.groupby(['Protein Accession #']).JDscore.idxmax()].sort_values(by=["JDscore"], ascending=False)
    
    #add PSM# count for each L_ID
    
#     psms_count=peptideDF.groupby(['Protein Accession #'])[['Scan#','PSM#','unique']].sum().reset_index()
#     psms_count2 =psms_count.rename(columns={"Scan#":"Total Peptide#",'unique':'Unique Peptide#'})
    
    psms_count2 = peptideDF.groupby('Protein Accession #').apply(my_agg)

    
    
    #remove peptide PSM# column later so that we can replace with protein PSM# column
    df_protein2 = df_protein.drop(['PSM#'], axis=1)
    df_protein_count = df_protein2.merge(psms_count2,how="outer",on="Protein Accession #")
    
    df_protein_count.rename(columns = {"Peptides":"Peptide of the Highest Score"},inplace=True)
    
    proteinFDR = calcFDR_df(df_protein_count)

    #drop Protein list column and just keep unique protein column
    df_protein_unique = df_protein_count.drop(['Protein'], axis=1)
    #map dictionary to Protein column
    
    
    df_protein_unique['Protein'] = df_protein_unique["Protein Accession #"].map(uniqProtAllProtDict)

#     For all protein split the Protein column 
    df_protein_countAll = tidy_split(df_protein_unique, "Protein", sep=',', keep=False)
    
    #drop Protein list column and just keep unique protein column
    df_protein_countAllFinal = df_protein_countAll.drop(["Protein Accession #"], axis=1)
    
    
    #rename the Protein column as Protein Accession #
    df_protein_countAllFinal.rename(columns={"Protein":"Protein Accession #"}, inplace=True)
    write_log("\nThe final protein FDR = {}\nThe unique protein level publication file is prepared".format(proteinFDR))
    
    #df_protein_count3 is for id_all_prot.txt and df_protein_count2 is id_uni_prot.txt
    return df_protein_unique,df_protein_countAllFinal
    

def ppmlFileReformat(ppmlFile):
    ppmlDF = fileToDF(ppmlFile)
    #we need to first make sure we sort the file by Protein_grp and Protein_sub_group
    
    ppmlDF["Protein_grp_num"] = ppmlDF["Protein_grp"].str.extract('(\d+)').astype(int)
    ppmlDF["subgroup"] = ppmlDF["Protein_sub_group"].astype("int")
    #ppmlDF.drop(["Peptides"], axis=1, inplace=True)
    #this is to make jump publication table like output
    ppmlDF.rename(columns={"Protein_grp":"group","PeptideSeqWithRealDelMass":"Peptide"}, inplace=True)

    #now sort by ascending order of first by Protein_grp_num and Protein_sub_group
    ppmlDFsorted = ppmlDF.sort_values(by=["Protein_grp_num","Protein_sub_group"], ascending=True)

    #make unique column 
    ppmlDFsorted["unique"] = np.where(ppmlDFsorted.Fate == "Unique","1","0")

    # #since JUMP -f cannot get all the genes of Tremble database (may be a bug). this part rescues GN of the same group and adds if missing
    # ppmlDFsortedGN = ppmlDFsorted.dropna(subset=["GN"])
    # grpToGeneDict = dict(zip(ppmlDFsortedGN.Protein_grp_num, ppmlDFsortedGN.GN))
    # ppmlDFsorted["GN"] = ppmlDFsorted.Protein_grp_num.map(grpToGeneDict)

    return ppmlDFsorted


def pitFileToRankDict(pitFile):
    pitDF = fileToDF(pitFile)
    rankDict = dict(zip(pitDF.ProteinName, pitDF.index))
    return rankDict


def id_prot_lookupFileParse(id_protFile):
    idProtDF = fileToDF(id_protFile)
    idProtDF["Peptide"] = idProtDF.L_ID.str.split(";",expand=True)[0]
    idProtDF_noDup = idProtDF.drop_duplicates()
    return idProtDF_noDup


def tidy_split(df, column, sep='|', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df



def rt_jdscore_filtering(df, rt_width, jdscore_width, delRT_cutoff, user_fdr=0.1):
    
    mean = np.mean(df.deltaRT_postcal)
    std = np.std(df.deltaRT_postcal)
    # global RT cutoff 2.5 SD of entire RT shift
    global_RT_cutoff_lt = mean - (std*delRT_cutoff)
    global_RT_cutoff_rt = mean + (std*delRT_cutoff)
    #binning based on JDscore
    #Big shift cutoff
    df_s_trim = df[df.deltaRT_postcal.between(global_RT_cutoff_lt,global_RT_cutoff_rt)]

    max_val_rt = np.round(np.max(df_s_trim.deltaRT_postcal),1)
    min_val_rt = np.round(np.min(df_s_trim.deltaRT_postcal),1)
    rt_bins = np.arange(min_val_rt-delRT_cutoff,max_val_rt+delRT_cutoff,delRT_cutoff)

    jdscore_bins = np.arange(0,1+jdscore_width,jdscore_width)

    #add RTbins and JDscore bins
    df_s_trim = bin_fixed(df_s_trim, rt_bins,"deltaRT_postcal")
    df_s_trim = bin_fixed(df_s_trim, jdscore_bins,"JDscore")
    
    jdscore_bins_sorted = sorted(set(df_s_trim.bins_JDscore))[::-1] #highest to lowest
    
    reqd_cols = df_s_trim.columns
    fdr_passed_bins = []
    fdr_failed_bins = []
    ctarget = 0
    cdecoy =0
    for jindex, jbin in enumerate(jdscore_bins_sorted):
        bin1_df = df_s_trim[df_s_trim.bins_JDscore == jbin]
    #     print(bin1_df)
        bins_deltaRT = sorted(set(bin1_df.bins_deltaRT_postcal)) #lowest to highest

            #for each bin compute FDR and see the number of targets
        target_track_list = []
        decoy_track_list = []
        bin_track_list = []
        # ctarget = 0
        # cdecoy =0
        for index, bins in enumerate(bins_deltaRT):
            # write_log ("Working for jdscore width of {} and rt width of {}\n".format(jbin, bins))
            boxdf = bin1_df[bin1_df.bins_deltaRT_postcal == bins]
            target = 0
            decoy = 0
            tar_decoy_dict = Counter(boxdf.Type)
#             print (tar_decoy_dict)
            if "Target" in tar_decoy_dict.keys():
                target = tar_decoy_dict["Target"]

            if "Decoy" in tar_decoy_dict.keys():
                decoy = tar_decoy_dict["Decoy"]

            if (target == 0) & (decoy == 0):
                print (" .. target = 0 and decoy = 0 for JDscore bin {} & deltaRT bin {}".format(jdscore_bins_sorted[0],bins))
                FDR=100
            elif target == 0:
                FDR = 100
            else:
                FDR = decoy/target*100
                if FDR > 100:
                    FDR = 100

            if FDR < user_fdr:
                target_track_list.append(target)
                decoy_track_list.append(decoy)
                bin_track_list.append(bins)
        
        if len(target_track_list) >=1:
            max_target = np.max(target_track_list)
            index_max_target = target_track_list.index(max_target)


            fdr_passed_index = [] # this will be updated for ascending and descending
            fdr_passed_rt_bins = []
            #compute cumulativeFDR descending index 
            for value in list(range(0,index_max_target+1))[::-1]:
                bin_target = target_track_list[value]
                bin_decoy = decoy_track_list[value]
                ctarget+=bin_target
                cdecoy+=bin_decoy

                cfdr = cdecoy/ctarget*100
                fdr_passed_index.append(value)
                fdr_passed_rt_bins.append(bin_track_list[value])
                if cfdr > user_fdr: # we are going 2 ways so we need to divide the user_fdr by 2
                    break   


            #compute cumulativeFDR ascending index 
            for value in list(range(index_max_target+1, len(target_track_list)-1)):
                bin_target = target_track_list[value]
                bin_decoy = decoy_track_list[value]
                ctarget+=bin_target
                cdecoy+=bin_decoy

                cfdr = cdecoy/ctarget*100
                fdr_passed_index.append(value)
                fdr_passed_rt_bins.append(bin_track_list[value])

                if cfdr > user_fdr: # we are going 2 ways so we need to divide the user_fdr by 2
                    break     

            bin1_pass = bin1_df[bin1_df.bins_deltaRT_postcal.isin(fdr_passed_rt_bins)]

            bin1_pass_FDR = FDR_Target_Decoy(bin1_pass,sortCol="JDscore")

            fdr_passed_bins.append(bin1_pass_FDR)
            bin1_fail = bin1_df[~bin1_df.bins_deltaRT_postcal.isin(fdr_passed_rt_bins)]
            fdr_failed_bins.append(bin1_fail)

    pass_bins_Df_FDR = pd.concat(fdr_passed_bins)
    fail_bins_Df = pd.concat(fdr_failed_bins)
    
    fail_bins_Df_FDR = FDR_Target_Decoy(fail_bins_Df,sortCol="JDscore")
    rescue_jdscore = fail_bins_Df_FDR.query("FDR<{}".format(user_fdr))
    
    final_psms_FDR = pd.concat([pass_bins_Df_FDR,rescue_jdscore])
    
    return final_psms_FDR