import scipy
from scipy import interpolate
import statsmodels.api as sm
import numpy as np  
import pandas as pd


def genLowessFunction(df, minv, maxv):
    X=list(df.ReferenceRT)
    Y=list(df.NewRT)
    
    x = [minv]+X+[maxv]
    y = [minv]+Y+[maxv]
    
    lowess = sm.nonparametric.lowess
    #Functions: lowess Fit a smooth nonparametric regression curve to a scatterplot
    Z = lowess(y,x,frac=0.05,it=3, return_sorted = False)
    #predicted y values using Z function
    y_ = interpolate.interp1d(x, Z, bounds_error=False)
    
    ### to generate the lowess curve use the xnew and ynew 
    #making a numpy array ranging from smallest reference/new library RT to largest
    #compute minimum and maximum value for the curve generation
#     minv = np.min(list(dfR_T.RT)+list(dfL_T.RT))
#     maxv = np.max(list(dfR_T.RT)+list(dfL_T.RT))
    xnew = np.arange(minv, maxv, 0.0001)
    
#     predicting ynew with the function y_
    ynew = y_(xnew)
    
    return y_, xnew, ynew

# def genLowessFunction(X, Y, minv, maxv, lowess_alpha):
#     ### to generate the lowess curve use the xnew and ynew 
#     #making a numpy array ranging from smallest reference/new library RT to largest
#     #compute minimum and maximum value for the curve generation
    
#     xnew = np.arange(minv, maxv, 0.001)
     
#     lowess = sm.nonparametric.lowess
#     #Functions: lowess Fit a smooth nonparametric regression curve to a scatterplot
# #     print (lowess_alpha)
#     Z = lowess(Y,X,frac=lowess_alpha,it=3, return_sorted = False)
#     #predicted y values using Z function
#     y_ = interpolate.interp1d(X, Z, bounds_error=False)
    

    
# #     predicting ynew with the function y_
#     ynew = y_(xnew)
    
#     return y_, xnew, ynew





def align_lowess(df, ref_run = "Yadav_B10"):
    # Initialization
    res = df.copy()
    colN = [c for c in df.columns if c.endswith("nPSMs")]
    colRt = df.columns.drop(colN)
    colRt = colRt.drop("key")

    all_cols = list(df.columns)
    refRunIdx = list(colRt).index(ref_run)

    # RT-alignment using a specific run as a reference (sequential alignment)
    # This approach is composed of three steps
    # 1. Set a specific run as a reference
    # 2. Align and calibrate one of remaining run against the reference using a LOESS model
    # 3. Update the reference by merging the old reference and the aligned/calibrated run (using weighted average)

    # Set the initial reference
    refIdx = refRunIdx
    ref = df[colRt[refIdx]].copy()
    refN = df[colN[refIdx]].copy()
    
    

    print("  {} is set to the initial reference for RT-alignment".format(colRt[refIdx]))

    # maximum value and minimum value 

    # Alignment and calibration of RT
    for i in range(len(colRt)):
        if i == refIdx:  # Skip for the reference run
            continue

            #get the adjacent run to reference for the lowess curve generation (get extreme points)
        #next run to ref
        next_ref = df[colRt[i]]
        
        # print (ref)

        minref = np.min(ref)
        maxref = np.max(ref)

        minnew = np.min(next_ref)
        maxnew = np.max(next_ref)


        minv = np.min([minref,minnew])
        maxv = np.max([maxref,maxnew])
        
    #     print (minv, maxv)
        print("  The RTs in {} are being aligned and calibrated".format(colRt[i]))

        # Prepare the modeling: find the shared peptides between the reference and the current run
        idx = (~ref.isna()) & (~next_ref.isna())
        ref_ = ref[idx]
        new_ = df[idx][colRt[i]]

        #optimize here
        max_refined_val, rt_cutoff, frac_val = lowess_optimization(ref_,new_, minv, maxv)

        print ("The optimized rt cutoff = {} and frac value = {}".format(rt_cutoff, frac_val))
#         #selecting best precursor for generating curve (tolerance within 2 minutes)
        idx_ = abs(ref_ - new_) < rt_cutoff

        x = ref_[idx_]
        y = new_[idx_]

#         print (len(x), len(y))
        # Build a LOESS model and calibrate RTs

        #create new x with added minref
        X = pd.concat([pd.Series([minv]), x, pd.Series([maxv])]) 
        Y = pd.concat([pd.Series([minv]), y, pd.Series([maxv])]) 
        

        mod, xpred, ypred = genLowessFunction(X, Y, minv, maxv,frac_val)

        res[colRt[i]] = mod(df[colRt[i]])

        # Update the reference by merging (taking the union of) the current reference and the calibrated run
        # There are three types of peptides (i.e., row indices)
        # 1. Peptides shared between the current reference and the calibrated run -> update the reference by taking the weighted average
        fullRt = pd.concat([ref[idx], res[colRt[i]][idx]], axis=1)
        fullN = pd.concat([refN[idx], res[colN[i]][idx]], axis=1)
        rt = (fullRt * fullN.values).sum(axis=1) / fullN.sum(axis=1)  # Reference as the weighted mean of RTs
        ref[idx] = rt
        refN[idx] = refN[idx] + res[colN[i]][idx]
        # 2. Peptides only found in the calibrated run -> replace the reference with the calibrated RTs
        idx = (ref.isna()) & (~res[colRt[i]].isna())
        ref[idx] = res[colRt[i]][idx]
        refN[idx] = res[colN[i]][idx]
        # 3. Peptides only found  in the current reference -> no action

    # Organization of the output dataframe
    # Calculation of the weighted standard deviation of RTs (https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf)
    M = (~res[colN].isna()).sum(axis=1)
    den = ((M - 1) / M) * res[colN].sum(axis=1)
    num = ((res[colRt].sub(ref, axis=0) ** 2) * res[colN].values).sum(axis=1)
    sdRt = np.sqrt(num / den)
    sdRt[den == 0] = 0
    res["SdRT"] = sdRt
    # # Calculation of the weighted average RTs
    # Calculation of the weighted average RTs
    res["AvgRT"] = ref  # In fact, the final reference is equivalent to the weighted average of the aligned/calibrated RTs
    
    return res,xpred, ypred, mod


############### Iterative appraoch ###############

def lowess_optimization(ref,new, minv, maxv):

    max_refined_val = 0
    rt_cutoff = 0
    frac = 0

    #generate lowess parameters here iteratively
    for rt_cutoff_val in range(1,5):
        for frac_val in  np.arange(0.1,1, 0.1):
            idx = abs(ref - new) < rt_cutoff_val
            x = ref[idx]
            y = new[idx]

            X = pd.concat([pd.Series([minv]), x, pd.Series([maxv])]) 
            Y = pd.concat([pd.Series([minv]), y, pd.Series([maxv])]) 

            mod, xpred, ypred = genLowessFunction(X, Y, minv, maxv,frac_val)
            cal_new = mod(new)

            idx_cal = abs(ref-cal_new) < 0.5
            refinedNo = len(ref[idx_cal])/len(ref)*100
            if refinedNo > max_refined_val:
                max_refined_val = refinedNo
                rt_cutoff = rt_cutoff_val
                frac = frac_val
    return max_refined_val, rt_cutoff, frac

############### END ##################
