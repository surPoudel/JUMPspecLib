#import sys
import pandas as pd
#import os
import numpy as np
#from tmtCorrection_mzXML import *
#import pandas as pd
import pyteomics
from pyteomics import mzxml,mass


# tmtReport = "TMT16" # this is parameter .. you can keep  
# mzFILE = "/Users/spoudel1/Desktop/TMT16_massCalibration/FTLD_Batch2_F50.mzXML"
# tol_max = 15 #this can be a parameter

# tmt_tag_mass = {"TMT0" :224.1524779, "TMT2" :225.1558327,
                # "TMT6" :229.1629321, "TMT7" :229.1629321,
               # "TMT8" :229.1629321, "TMT9" :229.1629321,
                # "TMT10" :229.1629321, "TMT11" :229.1629321,
               # "TMT16":304.2071453}



#this is the main function
#usage: firstSearchCorrection(mzFILE,tmt, tol_max=15 (optional)) returns you 2 dictionaries
# the input is mzXML file and tmt reporter ions list

def firstSearchCorrection(mzFILE,tmt, tol_max):
    
    tmt_tag_mass = {"TMT0" :224.1524779, "TMT2" :225.1558327,
                    "TMT6" :229.1629321, "TMT7" :229.1629321,
                   "TMT8" :229.1629321, "TMT9" :229.1629321,
                    "TMT10" :229.1629321, "TMT11" :229.1629321,
                   "TMT16":304.2071453}
    if len(tmt)==1:
        tmtReport = "TMT0"
    else:
        tmtReport = "TMT%d" % len(tmt)
    
    proton = 1.00727646677 #proton monoisotopic mass
    y1_Arg = str(mass.calculate_mass("R")+(1*proton)) #y1 ion with Arginine 
    y1_tmt_Lys = str(mass.calculate_mass("K")+(1*proton)+tmt_tag_mass[tmtReport]) #y1 ion with Lysine and TMT ion
    
    np_arr1,np_arr2,reportAllShift,lysineAllShift,arginineAllShift,massErrors = all_scans_TMTcorrection(mzFILE,tmt, y1_tmt_Lys, y1_Arg, tol_max)
    # superListMassShift = reportAllShift+lysineAllShift+arginineAllShift
    # correctionFactor = np.mean(superListMassShift)
    pos = np.nonzero(massErrors<=tol_max)[0]# get the mass shift with valid values (i.e. ms scans have references ions)
    if len(pos)==0:
        correctionFactor = 0
    else:
        correctionFactor = np.mean(massErrors[pos])
    
    run = mzFILE.split("/")[-1].split(".")[0]
    phage3plot = {}
    
    phage3plotPrecursoMass = {}
    
    for val in np_arr2:
        scan = str(val[0])
        exp_mz_list = val[2].tolist()
        intensity_list = val[3].tolist()
        precMZ = val[5][0]["precursorMz"]
        precIntensity = val[5][0]["precursorIntensity"]
        spectrum = run+"."+scan
        
        
        exp_mz_listC = massCorrectionFunction(exp_mz_list, massError=float(correctionFactor))
        precMZCorr = massCorrectionFunction([precMZ], massError=float(correctionFactor))

        valuePlotParent = [exp_mz_list,exp_mz_listC, intensity_list]
        phage3plot[spectrum] = valuePlotParent
        phage3plotPrecursoMass[spectrum] = [precMZ,precMZCorr[0],precIntensity]
        
    return phage3plot,phage3plotPrecursoMass

def all_scans_TMTcorrection(mzFILE,tmt, y1_tmt_Lys, y1_Arg,tol_max):
    np_arr1,np_arr2 = mzFileToNumpyArr(mzFILE) #getting the dataframe in numpy array form
#     print (np_arr2)
    correctionFactor_tmt = []
    correctionFactor_y1Lys_tmt = []
    correctionFactor_y1Arg = []
    
    # if the ms scan does not have the reference ions, the mass shift will be 2*tol_max
    # mass shift for each type (tmt, y1Lys_tmt, y1Arg)
    massErrors_tmt = np.zeros((len(np_arr2),1))+2*tol_max
    massErrors_y1Lys_tmt = np.zeros((len(np_arr2),1))+2*tol_max
    massErrors_y1Arg = np.zeros((len(np_arr2),1))+2*tol_max
    i=0
    for val in np_arr2:
        scan = str(val[0])
        exp_mz_list = val[2].tolist()
        intensity_list = val[3].tolist()
        #correctionFactor = massShiftCalculator(exp_mz_list,intensity_list,tmt,tol_max)# all tmt ions are used
        correctionFactor = massShiftCalculator(exp_mz_list,intensity_list,[tmt[0]],tol_max)# only tmt 126 is used
        correctionFactor_tmt+=correctionFactor
        correctionFactor_y1tmtLys = massShiftCalculator(exp_mz_list,intensity_list,[y1_tmt_Lys],tol_max)
        correctionFactor_y1Lys_tmt+=correctionFactor_y1tmtLys
        correctionFactor_y1R = massShiftCalculator(exp_mz_list,intensity_list,[y1_Arg],tol_max)
        correctionFactor_y1Arg+=correctionFactor_y1R
        
        # mass shift for each type (tmt, y1Lys_tmt, y1Arg)
        if len(correctionFactor)>0:
            massErrors_tmt[i]=np.mean(correctionFactor)
        if len(correctionFactor_y1tmtLys)>0:
            massErrors_y1Lys_tmt[i]=np.mean(correctionFactor_y1tmtLys)
        if len(correctionFactor_y1R)>0:
            massErrors_y1Arg[i]=np.mean(correctionFactor_y1R)
        i=i+1
    
    # select the type with the largest number of ms scans
    nvalidppm = [len(np.nonzero(massErrors_tmt<=tol_max)[0]),len(np.nonzero(massErrors_y1Lys_tmt<=tol_max)[0]),len(np.nonzero(massErrors_y1Arg<=tol_max)[0])]
    ix=np.argmax(nvalidppm)
    if ix==0:# tmt
        massErrors=massErrors_tmt
    elif ix==1:# y1Lys_tmt
        massErrors=massErrors_y1Lys_tmt
    else:# y1Arg
        massErrors=massErrors_y1Arg
    
    return np_arr1,np_arr2,correctionFactor_tmt,correctionFactor_y1Lys_tmt,correctionFactor_y1Arg,massErrors

def mzFileToNumpyArr(mzFILE):
    mzXcols = ['num', 'msLevel', 'm/z array', 'intensity array', 'retentionTime', 'precursorMz']
    
    x1 = pyteomics.mzxml.read(mzFILE)  #reading mzXML file using pyteomics
    dfMz = pd.DataFrame([x for x in x1])  #dataframe of the mzXML file
    #print (dfMz.columns)
    df = dfMz[mzXcols]
    ms1 = df.loc[df.msLevel==1]     #ms1 level scans
    np_arr1 = ms1.to_numpy()
    ms2 = df.loc[df.msLevel==2]     #ms2 level scans
    np_arr2 = ms2.to_numpy()
    return np_arr1,np_arr2

def massShiftCalculator(exp_mz_list,intensity_list,standardIonsList,tol_max=15): # it will be used for both MS1 and MS2, the highest peak will be selected
    checkMinReporter = standardIonsList[0]
    checkMaxReporter = standardIonsList[-1]
    
    min_exp_mzCheck = float(checkMinReporter)-(tol_max/100)
    max_exp_mzCheck = float(checkMaxReporter)+(tol_max/100)
    
    calibrationList = []
    calibrationListInt = []
    
    for i, val in enumerate(exp_mz_list):
        if min_exp_mzCheck <= float(val) <= max_exp_mzCheck:
            calibrationList.append(float(val))
            calibrationListInt.append(intensity_list[i])
    
#     print (calibrationList)
    reporterDict = {}            
    for reporters in standardIonsList:
        for index, masses in enumerate(calibrationList):
        
            massshift = ppmCalc(float(reporters), float(masses))
            if abs(massshift) > tol_max:
                pass
            else:
                if reporters not in reporterDict.keys():
                    reporterDict[reporters] = [[massshift,calibrationListInt[index],masses]]
                else:
                    reporterDict[reporters].append([massshift,calibrationListInt[index],masses])
    
    massShiftList = []
    if len(reporterDict)>0:
        for key in reporterDict.keys():
            intVal=0
            for subVal in reporterDict[key]:
                if subVal[1] >= intVal:
                    intVal = subVal[1]
                    massShiftKeep = subVal[0]
            massShiftList.append(massShiftKeep)
    
    return massShiftList

def massShiftCalculator_all(exp_mz_list,intensity_list,standardIonsList,tol_max=15): # it will be used for both MS1 and MS2, all peaks will be selected
    checkMaxReporter = standardIonsList[-1]
    checkMinReporter = standardIonsList[0]
    
    min_exp_mzCheck = float(checkMinReporter)-(tol_max/100)
    max_exp_mzCheck = float(checkMaxReporter)+(tol_max/100)
    
    
    calibrationList = []
    massShiftList = []
    
    for i, val in enumerate(exp_mz_list):
        if min_exp_mzCheck <= float(val) <= max_exp_mzCheck:
            calibrationList.append(float(val))
    
#     print (calibrationList)
    for reporters in standardIonsList:
        for index, masses in enumerate(calibrationList):
        
            massshift = ppmCalc(float(reporters), float(masses))
            if abs(massshift) > tol_max:
                pass
            else:
                massShiftList.append(massshift)
       
    return massShiftList

def calibratedMass(mz, massError):
    calibMass = mz/(1+(massError/1e6))
#     calibMass = mz - ((massError*mz)/1e6)
    return calibMass

def massCorrectionFunction(exp_list, massError):
    massCorrectedList = []
    for mz in exp_list:
        calibMass = calibratedMass(mz, massError)
        massCorrectedList.append(calibMass)
    return massCorrectedList

def ppmCalc(a, b, tol=10):
    a = float(a)
    b = float(b)
  #   massError = abs(a-b)*1e6/a
    massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
    return float(massError)

def MS2MassCorrection(mzFILE,tmt, tol_max):
    
    tmt_tag_mass = {"TMT0" :224.1524779, "TMT2" :225.1558327,
                    "TMT6" :229.1629321, "TMT7" :229.1629321,
                   "TMT8" :229.1629321, "TMT9" :229.1629321,
                    "TMT10" :229.1629321, "TMT11" :229.1629321,
                   "TMT16":304.2071453}
    if len(tmt)==1:
        tmtReport = "TMT0"
    else:
        tmtReport = "TMT%d" % len(tmt)
    
    proton = 1.00727646677 #proton monoisotopic mass
    y1_Arg = str(mass.calculate_mass("R")+(1*proton)) #y1 ion with Arginine 
    y1_tmt_Lys = str(mass.calculate_mass("K")+(1*proton)+tmt_tag_mass[tmtReport]) #y1 ion with Lysine and TMT ion
    
    np_arr1,np_arr2,reportAllShift,lysineAllShift,arginineAllShift,massErrors = all_scans_TMTcorrection(mzFILE,tmt, y1_tmt_Lys, y1_Arg, tol_max)
    MS1_index=np.zeros((len(np_arr1),5)) # MS1 scan, MS1 rt, MS1 peak num, baseline, massshift
    MS2_index=np.zeros((len(np_arr2),8)) # MS1 scan,MS1 rt,MS2 scan,m/z,z,Fragtype,MS2 peak num,massshift
    i=0
    for val in np_arr1:
        MS1_index[i,0]=int(val[0])
        i=i+1
    i=0
    for val in np_arr2:
        MS2_index[i,2]=int(val[0])
        i=i+1
    
    # superListMassShift = reportAllShift+lysineAllShift+arginineAllShift
    # correctionFactor = np.mean(superListMassShift)
    pos = np.nonzero(massErrors<=tol_max)[0]# get the mass shift with valid values (i.e. ms scans have references ions)
    if len(pos)==0:
        correctionFactor = 0
    else:
        correctionFactor = np.mean(massErrors[pos])
    # print([1,len(reportAllShift),np.mean(reportAllShift)])
    # print([2,len(lysineAllShift),np.mean(lysineAllShift)])
    # print([3,len(arginineAllShift),np.mean(arginineAllShift)])
    
    return MS1_index,MS2_index,massErrors

def MS1MassCorrection(mzFILE, referenceMasses, tol_max):
    
    np_arr1,np_arr2 = mzFileToNumpyArr(mzFILE) #getting the dataframe in numpy array form
    MS1_index=np.zeros((len(np_arr1),5)) # MS1 scan, MS1 rt, MS1 peak num, baseline, massshift
    MS2_index=np.zeros((len(np_arr2),8)) # MS1 scan,MS1 rt,MS2 scan,m/z,z,Fragtype,MS2 peak num,massshift
    i=0
    for val in np_arr1:
        MS1_index[i,0]=int(val[0])
        i=i+1
    i=0
    for val in np_arr2:
        MS2_index[i,2]=int(val[0])
        i=i+1
    
    massErrors = np.zeros((len(np_arr1),1))+2*tol_max
    i=0
    for val in np_arr1:
        scan = str(val[0])
        exp_mz_list = val[2].tolist()
        intensity_list = val[3].tolist()
        correctionFactor_ref = massShiftCalculator(exp_mz_list,intensity_list,referenceMasses,tol_max=15)
        if len(correctionFactor_ref)>0:
            massErrors[i]=np.mean(correctionFactor_ref)
        i=i+1
        
    return MS1_index,MS2_index,massErrors
