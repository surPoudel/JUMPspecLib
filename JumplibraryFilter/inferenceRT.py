import pandas as pd
import pyteomics
from pyteomics import mzxml
import numpy as np
from logFunctions import *
from scipy import interpolate
import statsmodels.api as sm

def rt_inference(mzFILE, idtxtdf):
    #parse mzxml file to dataframe
    #initializing the list that updates the infered rt
    final_rt_list = []
    #this function converts the mzxml file to dataframe and numpy array
    dfMz = mzXMLToNumpyArray(mzFILE)
    
    #grouping of scans based on peptide sequence
    peptideGroupedScanDF = peptideScanGrouping(idtxtdf)
    pepScanCols = list(peptideGroupedScanDF.columns)
    np_arr_pep = peptideGroupedScanDF.to_numpy()
    for pep_row in np_arr_pep:
        scanList = pep_row[pepScanCols.index("scan")]
        
        pepSp_dfMz = dfMz.loc[dfMz.num.astype("int").isin(scanList)]
        inferRT = inferPeptideRT(pepSp_dfMz)
        final_rt_list.append(inferRT)
    peptide_RT_dictionary = dict(zip(peptideGroupedScanDF.Peptide,final_rt_list))
#     peptideGroupedScanDF["peptideLevelRt"] = final_rt_list
    
    return peptide_RT_dictionary

#this function groups the scan using peptides as key
#this grouping is used later to compute weighted RT
def peptideScanGrouping(idtxtdf):
    #dfFirst = pd.read_csv(idtxt, delimiter="\t")
    scanList = idtxtdf.groupby(["Peptide"])["scan"].apply(list).reset_index()
    peptideGroupedScanDF = pd.DataFrame(scanList)
    return peptideGroupedScanDF

#this function converts mzxml file to dataframe using pyteomics
def mzXMLToNumpyArray(mzFILE):
    x1 = pyteomics.mzxml.read(mzFILE)  #reading mzXML file using pyteomics
    dfMz = pd.DataFrame([x for x in x1])  #dataframe of the mzXML file
    return dfMz

#((𝑅𝑇_1 𝐼_1+𝑅𝑇_2 𝐼_2+𝑅𝑇_3 𝐼_3+𝑅𝑇_4 𝐼_4 ))/(𝐼_1+𝐼_2+𝐼_3+𝐼_4 )
def weightedRT(rt_list, intensity_list):
    a = np.array(rt_list)
    b = np.array(intensity_list)
    product = a*b
    num = np.sum(product)
    den = np.sum(b)
    return num/den

#this function extract the survey precursor retention time and intensity and uses a function called weightedRT
#to calculate the weighted RT based on the extracted intensity and retension time for all the scans that are 
#group togenther by grouping functions
def inferPeptideRT(dfMz):
    intensity_List = []
    mz_cols = list(dfMz.columns)
    np_arr = dfMz.to_numpy()
    for row in np_arr:
        scan = row[mz_cols.index("num")]
        surveyRT = row[mz_cols.index("retentionTime")]
        precursorInfo = row[mz_cols.index("precursorMz")]
        precIntensity = precursorInfo[0]['precursorIntensity']
        if precIntensity == 0:
            precIntensity=row[mz_cols.index("basePeakIntensity")]
  
        intensity_List.append(precIntensity)
    intensity = np.array(intensity_List)
    
    rt_array = np.array(dfMz.retentionTime)
    if np.sum(intensity) == 0:
        print (scan, rt_array, intensity) 
    rt_weight = weightedRT( rt_array,intensity)
    return rt_weight


def RT_inference_MS2based(exp_mzxml, printDF2Rank1):
    #main program
    peptide_RT_dictionary = rt_inference(exp_mzxml, printDF2Rank1)

    printDF2Rank1["peptide_RT"] = printDF2Rank1["Peptide"].map(peptide_RT_dictionary)

    #delta RT is library RT - infered RT
    # printDF2Rank1["deltaRT"] = printDF2Rank1.RT.astype("float") - printDF2Rank1.peptide_RT.astype("float")
    # printDF2Rank1["abs_delRT"] = printDF2Rank1.deltaRT.abs()
    
    # printDF2Rank1.to_excel(outputFolder+"/Library_Search_Rank1.xlsx", index=None)

    # deltaRT_mean = np.mean(printDF2Rank1.deltaRT)
    # deltaRT_std = np.std(printDF2Rank1.deltaRT)
    
    # write_log("  The mean delta RT between the library and the fraction = {} min".format(np.round(deltaRT_mean,2)))
    # write_log("  The standard deviation delta RT between the library and the fraction = {} min".format(np.round(deltaRT_std,2)))
    
    
#########RT calibration using the reference RT#############
#We use Lowess for this purpose
# RT-alignment using a specific run as a reference (sequential alignment)
# This approach is composed of three steps
# 1. Set a specific psms (HQ) that has Library RT in it
# 2. Use the Library RT (RT column) to make a LOWESS curve
# 3. Align and calibrate all the PSMS RT using Lowess curve


def genLowessFunction(df, minv, maxv):
    X=list(df.RT)
    Y=list(df.peptide_RT)
    
    lowess = sm.nonparametric.lowess
    #Functions: lowess Fit a smooth nonparametric regression curve to a scatterplot
    Z = lowess(Y,X,frac=0.05,it=3, return_sorted = False)
    #predicted y values using Z function
    y_ = interpolate.interp1d(X, Z, bounds_error=False)
    
    ### to generate the lowess curve use the xnew and ynew 
    #making a numpy array ranging from smallest reference/new library RT to largest
    #compute minimum and maximum value for the curve generation
#     minv = np.min(list(dfR_T.RT)+list(dfL_T.RT))
#     maxv = np.max(list(dfR_T.RT)+list(dfL_T.RT))
    xnew = np.arange(minv, maxv, 0.0001)
    
#     predicting ynew with the function y_
    ynew = y_(xnew)
    
    return y_, xnew, ynew



def calibrateRT(searchFile, mzxml):
    df_s = pd.read_csv(searchFile)
    #this is key to add RT value to each fractions
    df_s["key_psms"] = df_s.L_ID+"."+df_s["Outfile"]
    df_s_all = df_s.copy()

    #Decide on the set of good quality precursors to infer RT.
    #Select all targets that has JDscore > 0.8
    df_s = df_s.loc[(df_s.Type == "Target") & (df_s.JDscore > 0.8)]
    #Target and decoy High score

    #df_s = df_s.loc[(df_s.Type == "Target") & (df_s.JDscore > 0.9)]
    #df_s = df_s.loc[df_s.JDscore > 0.9]


    write_log("Inferring RT for each PSMS in the fraction {}\n".format(mzxml))
    RT_inference_MS2based(mzxml, df_s)


    #compute minimum and maximum value for the curve generation
    minv = np.min(list(df_s.peptide_RT)+list(df_s.RT))
    maxv = np.max(list(df_s.peptide_RT)+list(df_s.RT))


    print ("\n\n****** RT RANGE INFORMATION PRIOR TO CALIBRATION *********\n")

    print ("Minimum RT value (min) for Library = ",np.min(df_s.RT))
    print ("Maximum RT value (min) for Library = ",np.max(df_s.RT))

    print ("Minimum inferred RT value (min) for all psms = ",np.min(df_s.peptide_RT))
    print ("Maximum inferred RT value (min) for all psms = ",np.max(df_s.peptide_RT))

    print ("\n\n*************** LOWESS CURVE INFORMATION ******************\n")

    print ("\nFor lowess curve generation, RT ranging from ", minv, " and ", maxv, " minutes is considered")

    #generate lowess fucntion (f) using overlapped peptides modelRTDF
    #xnew and ynew are required for the lowess curve plot later in the scatterplot
    print ("Generating lowess curve to allign New library RT with Reference RT\n")

    f,xnew,ynew = genLowessFunction(df_s, minv, maxv)

    write_log ("Calibrating the New Library RT using the LOWESS curve\n")

    #ALl the new RTs are calibrated by applying the fucntion
    df_s_all["peptide_RT"] = df_s_all.apply(lambda x: f(x.RT), axis=1)
    df_s_all.peptide_RT = df_s_all.peptide_RT.astype("float")

    df_s_all["deltaRT"] = df_s_all.RT.astype("float") - df_s_all.peptide_RT
    df_s_all["abs_delRT"] = df_s_all.deltaRT.abs()

    
    return df_s_all