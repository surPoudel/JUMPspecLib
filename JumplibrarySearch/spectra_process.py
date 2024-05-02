import numpy as np
# import pprint

import os
# import sys
import pickle
import pandas as pd
import re
# import time
import math
# from numba import jit


#### 0. def
# Change_ext
def Change_ext(original_fullfile,new_ext):
    # change the extension of fullfile
    [tpath,tname] = os.path.split(original_fullfile)
    [fname,ext] = os.path.splitext(tname)
    new_name = '%s%s' % (fname,new_ext)
    new_fullfile = os.path.join(tpath,new_name)
    return new_fullfile

# @jit(nopython=True)
# preprocess_topk_per_100mz
def preprocess_topk_per_100mz(exp_mz, exp_int, top_ions, binsize):# top_ions = 10 # take highest top 10 ions using intensity to rank
    # 1. cut off between (cutoffLowMass, cutoffHighMass)
    norm_1e2 = 1e2 # RT
    norm_1e4 = 1e4 # prec_MZ
    cutoffLowMass = 176 # Low Mass cutoff
    cutoffHighMass = 2000 # high mass limit
    exp_mz = np.array(exp_mz)
    exp_int = np.array(exp_int)
    ix = np.nonzero((exp_mz>=cutoffLowMass) & (exp_mz<=cutoffHighMass))[0] # cut off by range of cutoffLowMass and cutoffHighMass
    mz_list = exp_mz[ix].tolist()
    intensity_list = exp_int[ix].tolist()
    
    # 2. binning_mz_100
    top_mz_list = []
    top_intensity_list = []
    
    maxv = max(mz_list)+binsize
    minv = min(mz_list)-binsize
    bins = np.arange(minv, maxv, binsize) # fixed bin width 100 mz
    
    start_val = 0
    for x in range(1,len(bins)):
        sub_list_mz = [mzval for mzval in mz_list if bins[x-1] <= mzval < bins[x]] # subset mz list w
        if len(sub_list_mz) >=1: # some bins may not have any ions so we need to check the bin's element
            index_int = mz_list.index(sub_list_mz[-1]) # index for last value of sub list mz
            sub_list_int = intensity_list[start_val:index_int+1] # extract intensity until last value of sub list using intensity list (start and end index)
            start_val = index_int+1 # update start value for next round
            
            ind = np.argsort([-i for i in sub_list_int]) # sorting by descending order intensity
            r1 = np.argsort(ind)
            ix=np.nonzero(r1<top_ions)[0]
            
            # ind = np.argsort(sub_list_int) # sorting by descending order intensity
            # r1 = np.argsort(ind)
            # ix=np.nonzero(r1>=len(sub_list_int)-top_ions)[0]
            
            top_ion_lib=np.array(sub_list_mz)[ix].tolist()
            top_int_lib=np.array(sub_list_int)[ix].tolist()
            
            top_mz_list+=top_ion_lib  # append top 6 ions in the bin
            top_intensity_list+=top_int_lib # append top 6 ions intensity in the bin
    
    # 3. the intensity used are throughout normalized for the dot product calculation
    # simplifiedIonsDict = {"mz":top_mz_list,"intensity":top_intensity_list} # this dictionary stores mz and intensity after selecting top ions (example 10) after binning into 100 mz window
    mz=norm_1e4*np.array(top_mz_list)
    mz=[int(x)*1.0/norm_1e4 for x in mz]
    inten=np.array(top_intensity_list)
    
    # norm by the highest peak or the 2-fold of sum of peaks
    topinten_highest=max(inten)
    topinten_2sum=2*sum(inten)
    
    # norm by the highest peak
    inten=norm_1e2*inten/topinten_highest
    inten=[int(x) for x in inten]
    
    # save norm_factor=topinten_highest/topinten_2sum
    norm_factor=int(norm_1e4*topinten_highest/topinten_2sum)
    
    return mz,inten,norm_factor

# get_spec_df_from_pkl
def get_spec_df_from_pkl(specLib_pkl, top_ions, binsize):
    
    norm_1e2 = 1e2 # RT
    norm_1e4 = 1e4 # prec_MZ
    
    # specLib_norm_pkl
    specLib_norm_pkl = Change_ext(specLib_pkl,'_top'+str(int(top_ions))+'_per_100mz.pkl')
    if os.path.isfile(specLib_norm_pkl)==True:
        file = open(specLib_norm_pkl, "rb")
        expDF = pickle.load(file)
        file.close()
        return expDF
    
    # specLib_pkl
    if os.path.isfile(specLib_pkl)==False:
        data = {"scan": []}
        expDF = pd.DataFrame(data)
        print('file does not exist: {}'.format(specLib_pkl))
        return expDF
    
    # expDF
    file = open(specLib_pkl, "rb")
    expDF = pickle.load(file)
    file.close()
    
    # prec_MZ, RT, norm_factor
    expDF["prec_MZ"] = expDF.apply(lambda x: int(norm_1e4*x.prec_MZ)*1.0/norm_1e4, axis=1)
    expDF["RT"] = expDF.apply(lambda x: int(norm_1e2*x.RT)*1.0/norm_1e2, axis=1)
    expDF["norm_factor"] = expDF.apply(lambda x: int(x.RT), axis=1)
    
    # mz_cols1, np_arr1
    mz_cols1 = list(expDF.columns)
    np_arr1 = expDF.to_numpy()
    
    # m/z, intensity, norm_factor
    ms2_mz = []
    ms2_int = []
    ms2_normf = []
    for tno in range(0,len(np_arr1)):
        row1=np_arr1[tno]
        
        # load pkl (lib_data.columns: Index(['scan', 'charge', '[M+H]+', 'prec_MZ', 'L_ID', 'L_peptide', 'L_protein', 'RT', 'm/z', 'intensity', 'norm_factor']) )
        exp_mz_list = row1[mz_cols1.index("m/z")]
        exp_inten_list = row1[mz_cols1.index("intensity")]
        
        # picks up top 10 ions from each bin (100 mz)
        [exp_mz_list,exp_inten_list,norm_factor] = preprocess_topk_per_100mz(exp_mz_list,exp_inten_list,top_ions,binsize)
        
        # m/z, intensity
        ms2_mz.append(exp_mz_list)
        ms2_int.append(exp_inten_list)
        ms2_normf.append(norm_factor)
    
    # m/z, intensity, norm_factor
    dict1 = {"m/z":ms2_mz,"intensity":ms2_int,"norm_factor":ms2_normf}
    expDF1 = pd.DataFrame.from_dict(dict1)
    expDF["m/z"] = expDF1["m/z"]
    expDF["intensity"] = expDF1["intensity"]
    expDF["norm_factor"] = expDF1["norm_factor"]
    
    # sort expDF by prec_MZ
    expDF = expDF.sort_values(by='prec_MZ')
    
    expDF.to_pickle(specLib_norm_pkl)
    
    return expDF

# get_spec_df_from_ms2
def get_spec_df_from_ms2(exp_ms2, top_ions, binsize, sim_mass=0.0):
    
    norm_1e2 = 1e2 # RT
    norm_1e4 = 1e4 # prec_MZ
    
    if os.path.isfile(exp_ms2)==False:
        data = {"scan": []}
        expDF = pd.DataFrame(data)
        print('file does not exist: {}'.format(exp_ms2))
        return expDF
    
    g = open(exp_ms2,"r")
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
    ms2_normf = []

    mz_list = []

    for line in lines:
        if "S\t" in line:
            if len(mz_list) >= 1:
                [mz_list,int_list,norm_factor] = preprocess_topk_per_100mz(mz_list,int_list,top_ions,binsize)
                ms2_mz.append(mz_list)
                ms2_int.append(int_list)
                ms2_normf.append(norm_factor)
            mz_list = []
            int_list = []
    #         print ("spectrum found")
            temp_line = line.strip().split("\t")
            scan = int(temp_line[1])
            scan_list.append(scan)
            precursorIon = float(temp_line[-1])+sim_mass
            precursorIon_list.append(precursorIon)
        elif "Z\t" in line:
    #         print ("Charge found")
            temp_line = line.strip().split("\t")
            charge = int(temp_line[1])
            charge_list.append(charge)

            MH = float(temp_line[-1])
            MH_list.append(MH)
        elif "I\tRetTime" in line:# from *.ms2
    #         print ("RT found")
            temp_line = line.strip().split("\t")
            RT = float(temp_line[2])
            RT_list.append(RT)

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
        
        elif "L\tRT" in line:# from *.splib
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

    [mz_list,int_list,norm_factor] = preprocess_topk_per_100mz(mz_list,int_list,top_ions,binsize)
    ms2_mz.append(mz_list)
    ms2_int.append(int_list)
    ms2_normf.append(norm_factor)

    dict1 = {"scan":scan_list,"charge":charge_list,"[M+H]+":MH_list,
                "prec_MZ":precursorIon_list,"L_ID":L_ID_list,"L_peptide":L_peptide_list,
                "L_protein":L_Protein_list,"RT":RT_list,"m/z":ms2_mz,"intensity":ms2_int,"norm_factor":ms2_normf}

    dict2 = {}
    for key in dict1.keys():
        if len(dict1[key]) != 0:
            dict2[key] = dict1[key]

    expDF = pd.DataFrame.from_dict(dict2)
    
    # prec_MZ, RT
    expDF["prec_MZ"] = expDF.apply(lambda x: int(norm_1e4*x.prec_MZ)*1.0/norm_1e4, axis=1)
    expDF["RT"] = expDF.apply(lambda x: int(norm_1e2*x.RT)*1.0/norm_1e2, axis=1)
    
    return expDF

# get_spectra_from_df
def get_spectra_from_df(expDF,nLib):
    
    mz_cols1 = list(expDF.columns)
    np_arr1 = expDF.to_numpy()
    
    # # all_premz
    # all_premz=np.array([row[3] for row in np_arr1])
    
    # # r1: rank of all_premz by ascending order (the last is the highest)
    # ix = np.argsort(all_premz)
    # r1 = np.argsort(ix)
    
    # query_spectra, exp_spec_ID
    query_spectra = []
    exp_spec_ID = [[] for _ in range(len(np_arr1))]  # Initialize with empty lists
    
    for rno in range(0,len(np_arr1)):
        row1=np_arr1[rno]
        # tno=np.nonzero(r1==rno)[0][0]
        # row1=np_arr1[tno]
        
        # load pkl (lib_data.columns: Index(['scan', 'charge', '[M+H]+', 'prec_MZ', 'L_ID', 'L_peptide', 'L_protein', 'RT', 'm/z', 'intensity', 'norm_factor']) )
        exp_prec_mz = float( row1[mz_cols1.index("prec_MZ")] )
        exp_prech = int( row1[mz_cols1.index("charge")] )
        
        exp_mz = row1[mz_cols1.index("m/z")]
        exp_inten = row1[mz_cols1.index("intensity")]
        
        # exp_spec_ID
        if nLib==1:# library spec
            L_RT = float( row1[mz_cols1.index("RT")] )
            L_ID = row1[mz_cols1.index("L_ID")]
            L_peptide = row1[mz_cols1.index("L_peptide")]
            L_protein = row1[mz_cols1.index("L_protein")]
            exp_spec_ID[rno] = [exp_prec_mz,exp_prech,L_RT,L_ID,L_peptide,L_protein]# exp_prec_mz,exp_prech,L_RT,L_ID,L_peptide,L_protein
        else:# query spec
            cur_scan = int(row1[mz_cols1.index("scan")])
            exp_spec_ID[rno] = [exp_prec_mz,exp_prech,cur_scan]# exp_prec_mz,exp_prech,cur_scan
        
        # cur_spec
        cur_spec = {}
        cur_spec["id"] = rno
        cur_spec["precursor_mz"] = exp_prec_mz
        cur_spec["peaks"] = []
        for idx, cur_mz in enumerate(exp_mz):
            cur_spec["peaks"].append([cur_mz, exp_inten[idx]])
        
        query_spectra.append(cur_spec)
    # pprint.pprint(query_spectra)
    
    return query_spectra,exp_spec_ID

# create_fragment_index
def create_fragment_index(expDF,specLib_pkl=''):
    
    if os.path.isfile(specLib_pkl)==True:
        # frag_index_file
        frag_index_file = Change_ext(specLib_pkl,'_frag_index.dat')
        if os.path.isfile(frag_index_file)==True:
            frag_index = pickle.load(open(frag_index_file,'rb'))
            return frag_index
    
    mz_cols1 = list(expDF.columns)
    np_arr1 = expDF.to_numpy()
    
    # frag_index
    frag_index = {}
    
    for rno in range(0,len(np_arr1)):
        row1=np_arr1[rno]
        
        # load pkl (lib_data.columns: Index(['scan', 'charge', '[M+H]+', 'prec_MZ', 'L_ID', 'L_peptide', 'L_protein', 'RT', 'm/z', 'intensity', 'norm_factor']) )
        exp_prec_mz = row1[mz_cols1.index("prec_MZ")]
        exp_mz = row1[mz_cols1.index("m/z")]
        
        # frag_index
        for ino in range(0,len(exp_mz)):
            cur_frag_mz = exp_mz[ino]
            cur_prec_mz = exp_prec_mz
            if cur_frag_mz not in frag_index.keys():
                frag_index.update({cur_frag_mz:{cur_prec_mz:[rno]}})
            else:
                if cur_prec_mz not in frag_index[cur_frag_mz].keys():
                    frag_index[cur_frag_mz][cur_prec_mz]=[rno]
                else:
                    frag_index[cur_frag_mz][cur_prec_mz].append(rno)
    
    if os.path.isfile(specLib_pkl)==True:
        pickle.dump(frag_index,open(frag_index_file,'wb'))
    
    return frag_index

# get_lib_candidates
def get_lib_candidates(frag_index,all_frag,exp_prec_mz,exp_mz,ms1_tol,ms2_tol,match_threshold):
    
    # norm_1e4 = 1e4 # prec_MZ
    candidates = {}
    candidate_list = []
    
    delta_mz1 = ms1_tol*exp_prec_mz*1e-6
    left_lim1=exp_prec_mz-delta_mz1#int(norm_1e4*(exp_prec_mz-delta_mz1))*1.0/norm_1e4
    right_lim1=exp_prec_mz+delta_mz1#int(norm_1e4*(exp_prec_mz+delta_mz1))*1.0/norm_1e4
    
    for ino in range(0,len(exp_mz)):
        delta_mz2 = ms2_tol*exp_mz[ino]*1e-6
        left_lim2 = exp_mz[ino]-delta_mz2#int(norm_1e4*(exp_mz[ino]-delta_mz2))*1.0/norm_1e4
        right_lim2 = exp_mz[ino]+delta_mz2#int(norm_1e4*(exp_mz[ino]+delta_mz2))*1.0/norm_1e4
        
        pos2=np.nonzero((all_frag>=left_lim2) & (all_frag<=right_lim2))[0]
        if len(pos2)==0:
            continue
        
        for cur_R_frag_mz in all_frag[pos2]:
            cur_prec = np.array([x for x in frag_index[cur_R_frag_mz].keys()])
            
            pos1=np.nonzero((cur_prec>=left_lim1) & (cur_prec<=right_lim1))[0]
            if len(pos1)==0:
                continue
            
            for cur_R_prec_mz in cur_prec[pos1]:
                for jno in range(0,len(frag_index[cur_R_frag_mz][cur_R_prec_mz])):
                    cur_L_ID=frag_index[cur_R_frag_mz][cur_R_prec_mz][jno]
                    if cur_L_ID not in candidates.keys():
                        candidates.update({cur_L_ID:[1]})
                    else:
                        candidates[cur_L_ID][0]+=1
    
    for cur_L_ID in candidates.keys():
        if candidates[cur_L_ID][0]>=match_threshold:
            candidate_list.append(cur_L_ID)
    
    return candidate_list

# entropy_sim_per_peak
def entropy_sim_per_peak(x1,x2):
    # entropy similarity
    if x1==0.0 or x2==0.0:
        y = 0.0
    else:
        y = (x1+x2)*math.log((x1+x2),2) - x1*math.log(x1,2) - x2*math.log(x2,2)
    return y

# get_similarity
def get_similarity(e_intens,r_intens):
    
    if np.sqrt( sum(e_intens*e_intens)*sum(r_intens*r_intens) )==0.0:
        e_sim=0.0
    else:
        e_sim = sum(e_intens*r_intens)/np.sqrt( sum(e_intens*e_intens)*sum(r_intens*r_intens) )
    
    return e_sim

def get_pep_seq_only(seq):
    # change the extension of fullfile
    ABC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    seq_only=''
    for i in range(0,len(seq)):
        if seq[i] in ABC:
            seq_only=seq_only+seq[i]
    
    return seq_only

def get_delta_cn_rank(psm_score):
    # get_delta_cn_rank
    nrow = np.size(psm_score,0)
    # ncol = np.size(psm_score,1)
    delta_cn=np.zeros((nrow,1))
    
    score_col = np.array([psm_score[ino][0] for ino in range(0,nrow)])
    idx = (-score_col).argsort() # descending order
    high2low_rank = np.argsort(idx)
    
    score_col.sort() # ascending order
    for ino in range(0,nrow):
        c1_score = psm_score[ino][0]
        ix = np.nonzero(score_col<c1_score)[0]
        if len(ix)==0:
            c2_score = 0.0
        else:
            c2_score = score_col[ix[-1]]
        
        if c1_score==0.0:
            delta_cn[ino][0] = 0.0
        else:
            delta_cn[ino][0] = (c1_score-c2_score)/c1_score
    
    return delta_cn,high2low_rank
