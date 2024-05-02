#!/usr/bin/python

import os
import sys
import fnmatch
import copy
import pickle
import math
from numpy import *
from statistics import *
from DeisotopeMS1Features import *

def Get_params(jump_params):
	# get params
	
	jumpParamsDict = storeJUMPParams(jump_params)
	
	params={}
	
	# [Preprocess]
	params['deisotoping_method'] = int(float(jumpParamsDict['deisotoping_method']))	# 1: brute force; 2: pattern matching via entire envelope; 3: step-wise via adjacent peaks; 4: debug
	
	# parameters used by both brute force and top 1 for TMT
	params['precursor_ion_considered'] = int(float(jumpParamsDict['precursor_ion_considered']))	# N (1-4): go back to MS1, start from the strongest peak of top-N clusters in the isolation window
	params['charge_considered'] = get_int_list(jumpParamsDict['charge_considered'])	# a list of all charges to be considered
	params['monoisotopic_ion_shift'] = get_int_list(jumpParamsDict['monoisotopic_ion_shift'])	# a list of monoisotopic ion mass shift to be considered
	params['MS1_considered'] = get_int_list(jumpParamsDict['MS1_considered'])	# a list of MS1 scans to be considered, -1: MS1 immediate before the MS2; 1: MS1 immediate after the MS2; -2: the 2nd MS1 before the MS2
	
	# parameters used by top 1 for TMT
	params['number_TMT_tags'] = get_int_list(jumpParamsDict['number_TMT_tags'])	# a list of TMT tags to be considered
	params['isotopic_pattern_cutoff'] = float(jumpParamsDict['isotopic_pattern_cutoff'])	# percentage threshold (relative to the strongest ion) that a theoretical ion to be considered
	params['max_num_ppi'] = int(float(jumpParamsDict['max_num_ppi']))	# 0 = disable; 1-10 = max precursor ions selected for mixed MS2 search
	params['percentage_ppi'] = float(jumpParamsDict['percentage_ppi'])	# minimal percentage of precursor peak intensity (ppi) when max_num_ppi = 0
	params['isolation_window'] = float(jumpParamsDict['isolation_window'])	# +/- (isolation_window)/2 based on MS2 isolation window (e.g. 1.6 m/z)
	params['isolation_window_offset'] = float(jumpParamsDict['isolation_window_offset'])	# +/- isolation_window_offset based on MS2 isolation window offset (e.g. 0.25 m/z)
	params['isolation_window_variation'] = float(jumpParamsDict['isolation_window_variation'])	# +/- isolation_window_variation based on MS2 isolation window offset (e.g. 0.25 m/z)
	params['interscanppm'] = float(jumpParamsDict['interscanppm'])	# mass tolerance for interscan precursor identification
	params['intrascanppm'] = float(jumpParamsDict['intrascanppm'])	# mass tolerance for intrascan isotopic decharging
	params['delta_Dscore'] = float(jumpParamsDict['delta_Dscore'])	# threshold within which multiple precursors may be reported; delta_Dscore = (Dscore_max - Dscore_current) / Dscore_max
	params['TMT_data'] = int(float(jumpParamsDict['TMT_data']))	# 0 = disable; 1 = enable
	params['add_Nterm_peptide'] = float(jumpParamsDict['add_Nterm_peptide'])	# TMT modification or other amine labeling (229.1629321 for TMT6-11, 304.2071453 for TMT16)
	
	params['output_format'] = int(float(jumpParamsDict['output_format']))	# output file format, 1: ms2, 2: dtas, 3: mgf
	params['simple_process'] = int(float(jumpParamsDict['simple_process'])) # if 'deisotoping_method=1 or 2', 1: simple process with 2 steps; 2: complex process with 6 steps; if 'deisotoping_method=0', 1: charge +2; 2: charges +2, +3
	params['w_wo_calibration'] = int(float(jumpParamsDict['w_wo_calibration'])) # mass correction, 1: only with calibration, 2: only without calibration, else: both with and without calibration
	params['nprocessor'] = int(float(jumpParamsDict['nprocessor']))	# no of multiprocessing nodes
	params['parallel_method'] = int(float(jumpParamsDict['parallel_method'])) # parallel method, 1: multiprocessing package, 2: submitting LSF jobs
	
	tmt=get_tmt(jumpParamsDict)
	params['tmt'] = tmt	# TMT reporter ion mass table
	
	window = float(jumpParamsDict['isolation_window'])/2
	isolation_offset = float(jumpParamsDict['isolation_window_offset'])
	isolation_variation = float(jumpParamsDict['isolation_window_variation'])
	params['half_width_left'] = window-isolation_offset+isolation_variation # the left half width of the isolation window, unit: Th
	params['half_width_right'] = window+isolation_offset+isolation_variation # the right half width of the isolation window, unit: Th
	
	params['ms1_normalization'] = 0 # the function of MS1 normalization, 0: off; 1: on
	params['snr_limit'] = 1 # signal to noise ratio cutoff
	params['r2base_limit'] = 0.02 # signal to base peak ratio cutoff
	params['unitdiff'] = 1.00335 # the unit difference of isotopic peaks, close to mass difference of C13 and C12
	params['nextend'] = 1.1 # extention of isolation window, unit: Th
	params['tsim'] = 0.7 # similarity cutoff
	#params['nprocessor'] = 4 # no of multiprocessing nodes
	
	params['target_ms2scan'] = [9046]#18072#24338#	# specify the ms2 scan# to debug if deisotoping_method=4
	params['target_psm'] = [] # list of PSMs if deisotoping_method=4
	[tpath,tname] = os.path.split(jump_params)
	target_fullfile=os.path.join(tpath,'target_PSMs.xls')
	if True==os.path.isfile(target_fullfile):
		dat_fullfile = ReadTargetPSM(target_fullfile)
		psm = pickle.load(open(dat_fullfile,'rb')) # psm: scan_ch_mz_sq_mod_score_ac
		psm_scan = []
		for i in range(0,len(psm)):
			psm_scan.append(psm[i][0])
		params['target_ms2scan'] = psm_scan
		params['target_psm'] = psm
	
	return params

def storeJUMPParams(paramFile):
	# get jump params
	dict1 = {}
	with open(paramFile,"r") as f:
		for line in f:
			if (line.startswith("#") == False) and ("=" in line.strip()):
	
				te_line = line.strip().split("=")
				key = te_line[0].strip()

				if "#" in te_line[1]:
					value = te_line[1].split("#")[0].strip()
				else:
					value = te_line[1].strip()
				dict1[key]= value
	return dict1

def get_int_list(c_str):
	# int str ("," seperate) to int list 
	s_str = c_str.split(",")
	
	c_list = []
	for ino in range(0,len(s_str)):
		c_list.append(int(float(s_str[ino])))
	
	return c_list

def get_tmt(jumpParamsDict):
	# get tmt
	TMT_data = int(float(jumpParamsDict['TMT_data']))
	nTMT = 0
	if TMT_data==1:
		if abs(float(jumpParamsDict['add_Nterm_peptide'])-229.1629321)<0.0001:
			nTMT = 11
		elif abs(float(jumpParamsDict['add_Nterm_peptide'])-304.2071453)<0.0001:
			nTMT = 16
	
	if 16==nTMT:
		tmt = ['126.1277259','127.1247608','127.1310808','128.1281157','128.1344356','129.1314705','129.1377905',
		'130.1348253','130.1411453','131.1381802','131.1445001','132.1415350','132.1478550','133.1448899',
		'133.1512098','134.1482447']
	elif 11==nTMT:
		tmt = ['126.1277259','127.1247608','127.1310808','128.1281157','128.1344356','129.1314705','129.1377905',
		'130.1348253','130.1411453','131.1381802','131.1445001']
	else:
		tmt = []
	
	return tmt

def ReadTargetPSM(txt_fullfile):
	# read PSMs and save as dat
	# check if dat_fullfile exists
	dat_fullfile = Change_ext(txt_fullfile,'.dat')
	if True==os.path.isfile(dat_fullfile):
		psm_scan_ch_mz_sq_mod_score_ac = pickle.load(open(dat_fullfile,'rb'))
		return dat_fullfile
	
	# open
	file1 = open(txt_fullfile,'r')
	
	# read
	key1 = '	'
	psm_scan_ch_mz_sq_mod_score_ac = []
	c_str = file1.readline()
	pno = 0
	#print("PSMs:")
	while 1:
		# get a PSM
		c_str = file1.readline()
		if not c_str:
			break
		pno = pno+1
		#print(pno, end="\r")
		
		line = c_str.strip()
		domains = line.split(key1)
		#scan
		c_scan = int(float(domains[0]))
		#ch
		c_ch = int(float(domains[1]))
		#mz
		c_mz = float(domains[2])
		#sq (L2I)
		c_sq = domains[3]
		#mod
		c_mod = domains[4]
		#score
		c_score = float(domains[5])
		#ac
		c_ac = domains[-1]
		
		psm_scan_ch_mz_sq_mod_score_ac.append((c_scan,c_ch,c_mz,c_sq,c_mod,c_score,c_ac))
	
	# close
	#print("\n")
	file1.close()
	
	# sort and save
	#psm_scan_ch_mz_sq_mod_score_ac.sort(key=lambda x:x[5], reverse=True)
	psm_scan_ch_mz_sq_mod_score_ac.sort(key=lambda x:x[0])
	
	pickle.dump(psm_scan_ch_mz_sq_mod_score_ac,open(dat_fullfile,'wb'))
	return dat_fullfile

def write_mgf(file1,title,charge,rt,premz,p_mz,p_in):
	# write mgf
	file1.write('BEGIN IONS\n')
	file1.write('TITLE=%s\n' % title)
	file1.write('CHARGE=%d+\n' % charge)
	file1.write('RTINMINUTES=%.3f\n' % rt)
	file1.write('PEPMASS=%.7f\n' % premz)
	pnum = len(p_mz)
	for i in range(0,pnum):
		file1.write('%.7f %.1f\n' % (p_mz[i],p_in[i]))
	file1.write('END IONS\n')
	return

def write_dtas(file1,title,charge,premz,p_mz,p_in):
	# write dtas
	file1.write('%s %.7f %d\n' % (title,premz*charge-(charge-1)*1.007276,charge))
	pnum = len(p_mz)
	for i in range(0,pnum):
		file1.write('%.7f ' % (p_mz[i]))
	file1.write('\n')
	for i in range(0,pnum):
		file1.write('%.1f ' % (p_in[i]))
	file1.write('\n')
	return

def write_ms2(file1,scan,charge,rt,premz,pre_monointen,pre_monoppi,p_mz,p_in):
	# write ms2
	file1.write('S	%d	%d	%.7f\n' % (scan,scan,premz))
	file1.write('I	RetTime	%.3f\n' % rt)
	#file1.write('I	SumPrecursorIntensity	%.1f\n' % pre_monointen)
	#file1.write('I	PPI	%.3f\n' % pre_monoppi)
	file1.write('Z	%d	%.7f\n' % (charge,premz*charge-(charge-1)*1.007276))
	pnum = len(p_mz)
	for i in range(0,pnum):
		file1.write('%.7f\t%.1f\n' % (p_mz[i],p_in[i]))
	return

def OneFormat(SourcePath,extend):
	# get the full filename and only filename of the extended
	matched_fullfiles = []
	matched_files = []
	dir_list = os.listdir(SourcePath)
	for d in dir_list:
		if os.path.isfile(os.path.join(SourcePath,d)):
			if fnmatch.fnmatch(d.lower(),extend.lower()):
				matched_fullfiles.append(os.path.join(SourcePath,d))
				matched_files.append(d)
	return matched_fullfiles,matched_files

def Check_Path(output_path):
	# check and create path
	if False==os.path.isdir(output_path):
		os.mkdir(output_path)
		if False==os.path.isdir(output_path):
			print('can not create: ' % output_path)
	return

def Change_ext(original_fullfile,new_ext):
	# change the extension of fullfile
	[tpath,tname] = os.path.split(original_fullfile)
	[fname,ext] = os.path.splitext(tname)
	new_name = '%s%s' % (fname,new_ext)
	new_fullfile = os.path.join(tpath,new_name)
	return new_fullfile

def get_cur_MS2POS(num_MS2,cur_nbatch,total_nbatch):
	# get cur_MS2POS
	
	# cur_MS2POS
	nlen = num_MS2//total_nbatch
	st=(cur_nbatch-1)*nlen
	if cur_nbatch<total_nbatch:
		tm=cur_nbatch*nlen
	else:
		tm=num_MS2
	cur_MS2POS=arange(st,tm)
	
	return cur_MS2POS

def Load_tIPV(tIPV_fullfile):
	# load tIPV
	dat_tIPV_fullfile = Change_ext(tIPV_fullfile,'.dat')
	if True==os.path.isfile(dat_tIPV_fullfile):
		tIPV = pickle.load(open(dat_tIPV_fullfile,'rb'))
		return tIPV
	
	# read
	fid = open(tIPV_fullfile,"r")
	ftxt = fid.read()
	fid.close()
	ftxt_parts = ftxt.split()# tIPV with Inten: monoMH-1,monoMH,monoMH+1,monoMH+2,monoMH+3,monoMH+4,monoMH+5
	ncol = 7
	nIPV = int(len(ftxt_parts)/ncol)# length of tIPV
	
	# tIPV
	tIPV = []
	for i in range(0,nIPV):
		inten0 = float(ftxt_parts[i*ncol])
		inten1 = float(ftxt_parts[i*ncol+1])
		inten2 = float(ftxt_parts[i*ncol+2])
		inten3 = float(ftxt_parts[i*ncol+3])
		inten4 = float(ftxt_parts[i*ncol+4])
		inten5 = float(ftxt_parts[i*ncol+5])
		inten6 = float(ftxt_parts[i*ncol+6])
		tIPV.append((inten0,inten1,inten2,inten3,inten4,inten5,inten6))
	
	# save
	pickle.dump(tIPV,open(dat_tIPV_fullfile,'wb'))
	return tIPV

def Getall_tIPV(nTMT):
	# all_tIPV
	cur_path = os.getcwd()
	tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_noTMT_0tag.txt')
	all_tIPV_0=Load_tIPV(tIPV_fullfile)
	if 16==nTMT:
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT16_1tag.txt')
		all_tIPV_1=Load_tIPV(tIPV_fullfile)
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT16_2tag.txt')
		all_tIPV_2=Load_tIPV(tIPV_fullfile)
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT16_3tag.txt')
		all_tIPV_3=Load_tIPV(tIPV_fullfile)
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT16_4tag.txt')
		all_tIPV_4=Load_tIPV(tIPV_fullfile)
	elif 11==nTMT:
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT11_1tag.txt')
		all_tIPV_1=Load_tIPV(tIPV_fullfile)
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT11_2tag.txt')
		all_tIPV_2=Load_tIPV(tIPV_fullfile)
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT11_3tag.txt')
		all_tIPV_3=Load_tIPV(tIPV_fullfile)
		tIPV_fullfile = os.path.join(cur_path,'tIPV','tIPV_TMT11_4tag.txt')
		all_tIPV_4=Load_tIPV(tIPV_fullfile)
	if nTMT!=0:
		all_tIPV={"0": all_tIPV_0,"1": all_tIPV_1,"2": all_tIPV_2,"3": all_tIPV_3,"4": all_tIPV_4}
	else:
		all_tIPV={"0": all_tIPV_0}
	
	return all_tIPV

def Get_ms1_mode(MS1_peaks,index):
	# Get ms1_mode
	
	ms1_mode=1
	maxmz=[]
	for pos in range(0,4):
		mz0=MS1_peaks[index[pos]-1:index[pos+1]-1,0]
		maxmz.append(int(max(mz0)/100+0.5))
	
	if abs(maxmz[0]-maxmz[2])<=4 and abs(maxmz[1]-maxmz[3])<=4 and abs(maxmz[0]-maxmz[1])>4:
		ms1_mode=2
	
	return ms1_mode

def SelectMZ(mz,inten,startmz,nPeakWidthL,nPeakWidthR):
	# select peaks by mz on one MS1 scan
	
	sel_mz=[]
	sel_inten=[]
	left=startmz-nPeakWidthL
	right=startmz+nPeakWidthR
	pos=nonzero((mz>=left) & (mz<=right))[0]
	if len(pos)==0:
		return [sel_mz,sel_inten]
	sel_mz=mz[pos]
	sel_inten=inten[pos]
	return [sel_mz,sel_inten]

def SelectPeaks(MS1_index,MS1_peaks,index,params,ms1_mode,cur_premz,cur_ms1pos):
	# select peaks in the window on -/+ nscans MS1 scans
	
	nHalfWidthL=params['half_width_left']
	nHalfWidthR=params['half_width_right']
	nextend=params['nextend']
	snr_limit=params['snr_limit']
	r2base_limit=params['r2base_limit']
	ms1_list=params['MS1_considered'] # a list of MS1 scans to be considered, e.g. -2,-1,1 (when counting order, it is -1,0,1)
	ms1_normalization=params['ms1_normalization'] # the function of MS1 normalization, 0: off; 1: on
	
	num_MS1=MS1_index.shape[0]
	nPeakWidthL=nHalfWidthL+nextend
	nPeakWidthR=nHalfWidthR+nextend
	left_lim=cur_premz-nHalfWidthL
	right_lim=cur_premz+nHalfWidthR
	maxpeaknum=int(1e5) #the inti peak number of a MS2 spectrum
	
	# sel_pos
	sel_pos=[]
	for ino in ms1_list:
		if ino<0:# a list of MS1 scans to be considered, e.g. ms1_list -2,-1,1, while the order is -1,0,1
			ino=ino+1
		c_pos = cur_ms1pos+ino*ms1_mode
		if c_pos>=0 and c_pos<=num_MS1-1:
			sel_pos.append(c_pos)
	sel_pos=array(sel_pos)
	
	rt_len=len(sel_pos)
	sel_pnum=zeros((rt_len,1))
	sel_mz=zeros((rt_len,maxpeaknum))
	sel_inten=zeros((rt_len,maxpeaknum))
	
	for tno in range(0,rt_len):
		pos=sel_pos[tno]
		#print([num_MS1,len(index),pos])
		#print([index[pos],index[pos+1]])
		mz0=MS1_peaks[index[pos]-1:index[pos+1]-1,0]
		inten0=MS1_peaks[index[pos]-1:index[pos+1]-1,1]
		
		# select
		[mz,inten]=SelectMZ(mz0,inten0,cur_premz,nPeakWidthL,nPeakWidthR)
		baseline = MS1_index[pos,3]
		
		# cut the noise
		if size(inten)>0:
			i=nonzero((mz>left_lim)&( mz<right_lim))[0]
			if len(i)!=0:
				topin=max(inten[i])
			else:
				topin=0
			IX=nonzero(inten>=r2base_limit*topin)[0]
			mz=mz[IX]
			inten=inten[IX]
		
			IX=nonzero(inten>=snr_limit*baseline)[0]
			mz=mz[IX]
			inten=inten[IX]
		
		# save
		cur_len=len(mz)
		sel_pnum[tno]=cur_len
		if cur_len>0:
			sel_mz[tno,0:cur_len]=mz
			if ms1_normalization==1:
				c_max=max(inten)
				inten=inten/c_max*100.0
			sel_inten[tno,0:cur_len]=inten
	
	sno=int(max(sel_pnum))
	if sno<maxpeaknum:
		sel_mz=sel_mz[:,0:sno]
		sel_inten=sel_inten[:,0:sno]
	
	IX=nonzero(sel_pnum>0)[0]
	sel_pnum=sel_pnum[IX]
	sel_pos=sel_pos[IX]
	sel_mz=sel_mz[IX,:]
	sel_inten=sel_inten[IX,:]
	
	return [sel_pnum,sel_pos,sel_mz,sel_inten]

def CheckAdjMS1scans(params,ms1_mode,cur_premz,cur_ms1pos,sel_pnum,sel_pos,sel_mz,sel_inten,ncase,logfile):
	
	ptola=params['intrascanppm']
	ptole=params['interscanppm']
	unitdiff=params['unitdiff']
	Tsim=params['tsim']
	nrunning_mode=params['deisotoping_method']
	charge_list = params['charge_considered'] # a list of all charges to be considered
	ms1_normalization=params['ms1_normalization'] # the function of MS1 normalization, 0: off; 1: on
	
	# topn of merged MS1 scans
	[mz_m,inten_m] = GetRefMS1(sel_pnum,sel_mz,sel_inten,ptole)
	if ncase==1:
		[isonum,isomz,isointen,chg] = GetAllClustersOnSingleMS1(mz_m,inten_m,charge_list,unitdiff,ptola*2,1)
	else:
		[isonum,isomz,isointen,chg] = GetAllClustersOnSingleMS1(mz_m,inten_m,charge_list,unitdiff,ptola*2,2)
	inum=size(isomz,0)
	if 0==inum:
		if nrunning_mode==4:
			logfile.write('\n----------------CheckAdjMS1scans_log----------------\n')
			logfile.write('no clusters on merged MS1 ')
			if ncase==1:
				logfile.write('for charge +1\n')
			else:
				logfile.write('for charges +2, +3, etc\n')
		return [sel_pnum,sel_pos,sel_mz,sel_inten]
	
	#topn = 1
	[isonum,isomz,isointen,chg]=FilterByPPI(isonum,isomz,isointen,chg,0,cur_premz,ptola*2)
	c_len=int(isonum)
	mz_n = isomz[0,0:c_len]
	inten_n = isointen[0,0:c_len]
	
	# align cur MS1 to topn of merged MS1 scans
	[mz,inten] = GetOneMS1(sel_pnum,sel_pos,sel_mz,sel_inten,ms1_mode,cur_ms1pos,0)
	[mz_0,inten_0] = AlignMS1byMergedTopn(mz,inten,mz_n,inten_n,ptole*2)
	X=nonzero(inten_0>0)[0]
	if len(X)<=2:
		if nrunning_mode==4:
			logfile.write('\n----------------CheckAdjMS1scans_log----------------\n')
			logfile.write('<=2 peaks of the top1 cluster on MS1[-1]\n')
		return [sel_pnum,sel_pos,sel_mz,sel_inten]
	
	if nrunning_mode==4:
		logfile.write('\n----------------CheckAdjMS1scans_log----------------\n')
		logfile.write('mz_0\tinten_0\n')
		for i in range(0,len(mz_0)):
			logfile.write('%.4f\t%.0f\n' % (mz_0[i],inten_0[i]))
	
	# compare to cur MS1 (topn)
	e_sim_all=[]
	for ino in range(0,len(sel_pos)):
		c_len=int(sel_pnum[ino])
		mz = sel_mz[ino,0:c_len]
		inten = sel_inten[ino,0:c_len]
		[mz_i,inten_i] = AlignMS1byMergedTopn(mz,inten,mz_n,inten_n,ptole*2)
		[mlenR,mintenR]=match_inten(inten_i,inten_0)
		if sum(inten_i)==0:
			sel_pnum[ino]=0
			e_sim_all.append(0)
		else:
			e_sim = sum(inten_i*inten_0)/sqrt( sum(inten_i*inten_i)*sum(inten_0*inten_0) )
			e_sim_all.append(e_sim)
			if (ms1_normalization==1 and e_sim<Tsim) or ( ms1_normalization==0 and ((e_sim<Tsim and mlenR<0.6) or mintenR>5) ):
				sel_pnum[ino]=0
			if nrunning_mode==4:
				cur_no = int((sel_pos[ino]-cur_ms1pos)/ms1_mode)
				logfile.write('\n')
				logfile.write('mz_%d\tinten_%d\n' % (cur_no,cur_no))
				for i in range(0,len(inten_i)):
					logfile.write('%.4f\t%.0f\n' % (mz_i[i],inten_i[i]))
				logfile.write('++sim_%dvs0>=%.2f++\n' % (cur_no,Tsim))
				logfile.write('%.3f\n' % (e_sim))
	e_sim_all=array(e_sim_all)
	
	# check the 2nd max sim
	IX=nonzero(e_sim_all<1)[0]
	if len(IX)>0:
		sim_max_2nd=max(e_sim_all[IX])
		if sim_max_2nd>=0.85:
			for ino in range(0,len(sel_pos)):
				if e_sim_all[ino]<(1-0.1)*sim_max_2nd:
					sel_pnum[ino]=0
	
	IX=nonzero(sel_pnum>0)[0]
	sel_pnum=sel_pnum[IX]
	sel_pos=sel_pos[IX]
	sel_mz=sel_mz[IX,:]
	sel_inten=sel_inten[IX,:]
	
	if nrunning_mode==4:
		x=[]
		for ino in range(0,len(sel_pos)):
			cur_no = int((sel_pos[ino]-cur_ms1pos)/ms1_mode)
			x.append(cur_no)
		logfile.write('\nmerged ms1 scans:\n')
		for i in range(0,len(x)):
			logfile.write(' %.0f\n' % (x[i]))
		logfile.write('\n')
	
	return [sel_pnum,sel_pos,sel_mz,sel_inten]

def GetRefMS1(sel_pnum,sel_mz,sel_inten,ptol):
	# get reference MS1 peaks
	
	# merge and sort by mz
	peaks=[]
	for ino in range(0,len(sel_pnum)):
		c_len=int(sel_pnum[ino])
		for jno in range(0,c_len):
			peaks.append((sel_mz[ino,jno],sel_inten[ino,jno]))
	peaks.sort(key=lambda x:x[0])
	
	# mz0, inten0
	mz0=[]
	inten0=[]
	for ino in range(0,len(peaks)):
		mz0.append(peaks[ino][0])
		inten0.append(peaks[ino][1])
	mz0=array(mz0)
	inten0=array(inten0)
	
	# r1: rank of inten0 by ascending order (the last is the highest)
	ix = argsort(inten0)
	r1 = argsort(ix)
	
	# merge through the rank (the highest the first), and sort by mz
	peaks=[]
	procflag=zeros((len(mz0),1))
	for rno in range(len(mz0)-1,-1,-1):
		ino=nonzero(r1==rno)[0][0]
		if 1==procflag[ino]:
			continue
		
		c_ptol=ptol*mz0[ino]*1e-6
		left=mz0[ino]-c_ptol
		right=mz0[ino]+c_ptol
		pos=nonzero((mz0>=left) & (mz0<=right))[0]
		c_mz=mz0[pos]
		c_inten=inten0[pos]
		peaks.append((sum(c_mz*c_inten)/sum(c_inten),max(c_inten))) #### max of inten, or sum of inten: here choose max
		#peaks.append((sum(c_mz*c_inten)/sum(c_inten),sum(c_inten))) #### max of inten, or sum of inten: here choose max
		procflag[pos]=1
	peaks.sort(key=lambda x:x[0])
	
	# mz, inten
	mz=[]
	inten=[]
	for ino in range(0,len(peaks)):
		mz.append(peaks[ino][0])
		inten.append(peaks[ino][1])
	mz=array(mz)
	inten=array(inten)
	
	return [mz,inten]

def GetAllClustersOnSingleMS1(mz,inten,charge_list,unitdiff,ptol,ncase):
	# Get Clusters On the Reference MS1
	
	#set iso
	mz_len=len(mz)
	isonum=array(40*mz_len*[0])#zeros((mz_len,1))
	isomz=zeros((40*mz_len,60))
	isointen=zeros((40*mz_len,60))
	chg=array(40*mz_len*[0])
	inum=0
	
	#set peak flag :1 for searched 0 for not
	ref_chg=max(charge_list)
	procflag=zeros((mz_len,ref_chg))
	if ncase==1:
		chg_range=arange(1,1+1)#only 1
	else:
		chg_range=array([i for i in charge_list if i!=1])#start from 2
	#print([min(mz), max(mz)])
	
	## determin the candidate clusters: identical m/z intervals
	for pno in range(mz_len):
		for charge in chg_range:
			if 1==procflag[pno,charge-1]:
				continue
			
			#1 init			
			tp_isonum=0
			tp_isomz=[]
			tp_isomz.append(mz[pno])
			tp_isointen=[]
			tp_isointen.append(inten[pno])
			startmz=mz[pno]
			deltamz=unitdiff/charge
			
			#2 find a cluster: identical mz intervals
			while 1:
				#search the interval peak
				c_mz=startmz+deltamz
				c_ptol=ptol*c_mz*1e-6
				left=c_mz-c_ptol
				right=c_mz+c_ptol
				pos=nonzero((mz>=left) & ( mz<=right))[0]
				# find the interval peak
				if len(pos)!=0:
					#judge by intensity ratio
					tp_isonum=tp_isonum+1
					ix=argmax(inten[pos])
					if ( (charge>2 and tp_isonum>=3) or (charge==2 and tp_isonum>=2) ) and tp_isointen[tp_isonum-1]/inten[pos[ix]]<0.4:# cut into two clusters
						tp_isonum=tp_isonum-1
						break
					tp_isomz.append(mz[pos[ix]])
					tp_isointen.append(inten[pos[ix]])
					startmz=tp_isomz[tp_isonum]
				else:
					break
			
			#3 judge
			if 0==tp_isonum:
				continue
			
			#4 save
			isonum[inum]=tp_isonum+1
			isomz[inum,0:tp_isonum+1]=tp_isomz
			isointen[inum,0:tp_isonum+1]=tp_isointen
			chg[inum]=charge
			inum=inum+1
			
			#5 write the flag
			for tno in range(tp_isonum+1):
				pos=nonzero(tp_isomz[tno]==mz)[0]
				procflag[pos,charge-1]=1
	
	if 0==inum:
		mz_len=0
		isonum=array(40*mz_len*[0])#zeros((mz_len,1))
		isomz=zeros((40*mz_len,60))
		isointen=zeros((40*mz_len,60))
		chg=array(40*mz_len*[0])
		return [isonum,isomz,isointen,chg]
	
	if inum<40*mz_len:
		isonum=isonum[0:inum]
		isomz=isomz[0:inum,:]
		isointen=isointen[0:inum,:]
		chg=chg[0:inum]
	
	return [isonum,isomz,isointen,chg]

def FilterByPPI(isonum,isomz,isointen,chg,ncase,cur_premz,ptol):
	# filter by PPI
	
	inum=size(isomz,0)
	c_ptol=ptol*cur_premz*1e-6
	c_tol=0.5
	isoclusterint=array([0.0]*inum)
	flag_ppm=zeros((inum,1))
	flag_Da=zeros((inum,1))
	for ino in range(inum):
		isoclusterint[ino]=sum(isointen[ino,0:isonum[ino]])
		c_isomz=isomz[ino,0:isonum[ino]]
		IX=nonzero((cur_premz-c_ptol<=c_isomz) & (c_isomz<=cur_premz+c_ptol) )[0]
		if len(IX)>0:
			flag_ppm[ino]=1
		IIX=nonzero((cur_premz-c_tol<=c_isomz) & (c_isomz<=cur_premz+c_tol) )[0]
		if len(IIX)>0:
			flag_Da[ino]=1
	
	if ncase==0:
		topn=1
		if inum<=topn:
			return [isonum,isomz,isointen,chg]
		
		x1=nonzero(flag_ppm==1)[0]
		if len(x1)>0:
			ix = argsort(isoclusterint[x1])
			r1 = argsort(ix)
			xa=nonzero(r1>=len(x1)-topn)[0]
			isonum=isonum[x1[xa]]
			isomz=isomz[x1[xa],:]
			isointen=isointen[x1[xa],:]
			chg=chg[x1[xa]]
		else:
			ix = argsort(isoclusterint)
			r1 = argsort(ix)
			x=nonzero(r1>=inum-topn)[0]
			isonum=isonum[x]
			isomz=isomz[x,:]
			isointen=isointen[x,:]
			chg=chg[x]
	elif ncase==1:
		topn=1
		if inum<=topn:
			return [isonum,isomz,isointen,chg]
		
		x1=nonzero(flag_ppm==1)[0]
		if len(x1)>0:
			ix = argsort(isoclusterint[x1])
			r1 = argsort(ix)
			xa=nonzero(r1>=len(x1)-topn)[0]
			
			ix = argsort(isoclusterint)
			r1 = argsort(ix)
			x=nonzero(r1>=inum-topn)[0]
			
			if x[0]!=x1[xa[0]] and flag_Da[x[0]]==1 and isoclusterint[x[0]]/isoclusterint[x1[xa[0]]]>5:# much higher than the selected precursor
				isonum=isonum[x]
				isomz=isomz[x,:]
				isointen=isointen[x,:]
				chg=chg[x]
			else:
				isonum=isonum[x1[xa]]
				isomz=isomz[x1[xa],:]
				isointen=isointen[x1[xa],:]
				chg=chg[x1[xa]]
		else:
			ix = argsort(isoclusterint)
			r1 = argsort(ix)
			x=nonzero(r1>=inum-topn)[0]
			isonum=isonum[x]
			isomz=isomz[x,:]
			isointen=isointen[x,:]
			chg=chg[x]
	elif ncase>=2:
		topn=int(ncase)
		if topn%2==0:
			topns=[int(topn/2),int(topn/2)]
		else:
			topns=[int(topn/2)+1,int(topn/2)]
		if inum<=topn:
			return [isonum,isomz,isointen,chg]
		
		x1=nonzero(flag_ppm==1)[0]
		x2=nonzero(flag_ppm==0)[0]
		if len(x1)>0 and len(x2)>0:
			ix = argsort(isoclusterint[x1])
			r1 = argsort(ix)
			xa=nonzero(r1>=len(x1)-topns[0])[0]
			
			ix = argsort(isoclusterint[x2])
			r1 = argsort(ix)
			xb=nonzero(r1>=len(x2)-topns[1])[0]
			
			x=[]
			for ino in range(0,len(xa)):
				x.append(x1[xa[ino]])
			for ino in range(0,len(xb)):
				x.append(x2[xb[ino]])
			x=array(x)
			x=sort(x)
			
			isonum=isonum[x]
			isomz=isomz[x,:]
			isointen=isointen[x,:]
			chg=chg[x]
		else:# case a: len(x1)>0 and len(x2)==0, case b: len(x1)==0
			ix = argsort(isoclusterint)
			r1 = argsort(ix)
			x=nonzero(r1>=inum-topn)[0]
			isonum=isonum[x]
			isomz=isomz[x,:]
			isointen=isointen[x,:]
			chg=chg[x]
	
	return [isonum,isomz,isointen,chg]

def GetOneMS1(sel_pnum,sel_pos,sel_mz,sel_inten,ms1_mode,cur_ms1pos,nshift):
	
	mz=[]
	inten=[]
	c_pos = cur_ms1pos+nshift*ms1_mode
	ino=nonzero(sel_pos==c_pos)[0]
	if len(ino)!=0:
		c_len=int(sel_pnum[ino])
		c_mz = sel_mz[ino,0:c_len]
		c_inten = sel_inten[ino,0:c_len]
		c_mz = c_mz.flatten()
		c_inten = c_inten.flatten()
		for jno in range(0,c_len):
			mz.append(c_mz[jno])
			inten.append(c_inten[jno])
	
	mz=array(mz)
	inten=array(inten)
	
	return [mz,inten]

def AlignMS1byMergedTopn(mz,inten,mz_n,inten_n,ptol):
	
	mz_t=[]
	inten_t=[]
	
	for ino in range(0,len(mz_n)):
		c_ptol=ptol*mz_n[ino]*1e-6
		left=mz_n[ino]-c_ptol
		right=mz_n[ino]+c_ptol
		pos=nonzero((mz>=left) & (mz<=right))[0]
		if len(pos)==0:
			mz_t.append(mz_n[ino])
			inten_t.append(0.0)
		else:
			mz_t.append(sum(mz[pos]*inten[pos])/sum(inten[pos]))
			inten_t.append(max(inten[pos]))
	
	mz_t=array(mz_t)
	inten_t=array(inten_t)
	
	return [mz_t,inten_t]

def match_inten(inten_i,inten_0):
	
	ratio=[]
	for ino in range(0,len(inten_i)):
		if inten_i[ino]>0 and inten_0[ino]>0:
			ratio.append(inten_i[ino]/inten_0[ino])
	
	X=nonzero(inten_0>0)[0]
	nonzero_0=len(X)
	
	mlenR=len(ratio)/nonzero_0
	
	if len(ratio)>0:
		mintenR=median(ratio)
	else:
		mintenR=100
	
	return [mlenR,mintenR]

def FilterByISO(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,ptol):
	# filter by isotopic relationship
	
	[mz_len,c]=isomz.shape
	if mz_len==1:
		return [isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp]
	
	flag=ones((mz_len,1))
	for ino in range(0,mz_len-1):
		i_isonum=isonum[ino]
		i_isomz=isomz[ino,0:i_isonum]
		i_isointen=isointen[ino,0:i_isonum]
		i_chg=chg[ino]
		for jno in range(ino+1,mz_len):
			j_isonum=isonum[jno]
			j_isomz=isomz[jno,0:j_isonum]
			j_isointen=isointen[jno,0:j_isonum]
			j_chg=chg[jno]
			
			c_ptol=ptol*j_isomz[0]*1e-6
			ix=nonzero(abs(j_isomz-i_isomz[0])<=ptol*i_isomz[0]*1e-6)[0]
			iix=nonzero(abs(i_isomz-j_isomz[0])<=c_ptol)[0]
			if abs(i_isomz[0]-j_isomz[0])<=c_ptol:
				if i_chg==j_chg and i_isointen[0]<=j_isointen[0]:
					flag[ino]=0
				elif i_chg==j_chg and i_isointen[0]>j_isointen[0]:
					flag[jno]=0
				elif mod(j_chg,i_chg)==0:
					flag[ino]=0
				elif mod(i_chg,j_chg)==0:
					flag[jno]=0
			elif len(ix)>0 and j_chg>i_chg and mod(j_chg,i_chg)==0:
				flag[ino]=0
			elif len(iix)>0 and i_chg>j_chg and mod(i_chg,j_chg)==0:
				flag[jno]=0
	IX=nonzero(flag==1)[0]
	isonum = isonum[IX]
	isomz = isomz[IX,:]
	isointen = isointen[IX,:]
	chg = chg[IX]
	Dscore = Dscore[IX]
	delta_Dscore = delta_Dscore[IX]
	suminten = suminten[IX]
	isogrp = isogrp[IX]
	
	return [isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp]

def GetTopPPI(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,max_num_ppi,percentage_ppi):
	# filter by PPI
	
	grps=unique(isogrp)
	isoclusterint=array([0.0]*len(grps))
	for ino in range(0,len(grps)):
		ix=nonzero(isogrp==grps[ino])[0]
		isoclusterint[ino]=sum(suminten[ix])
	
	topn=max_num_ppi
	if topn==0:
		topn=1
	if topn==1 and max(isoclusterint)/sum(isoclusterint)<percentage_ppi/100.0:
		topn=2
	
	ix = argsort(isoclusterint)
	r1 = argsort(ix)
	x=nonzero(r1>=len(grps)-topn)[0]
	top_grps=grps[x]
	
	monomz = []
	monochg = []
	monointen = []
	monoppi = []
	for ino in range(0,len(top_grps)):
		ix=nonzero(isogrp==top_grps[ino])[0]
		for jno in range(0,len(ix)):
			c_ppi = suminten[ix[jno]]/sum(isoclusterint)
			if c_ppi<0.1:
				continue
			monomz.append(isomz[ix[jno],0])
			monochg.append(chg[ix[jno]])
			monointen.append(suminten[ix[jno]])
			monoppi.append(c_ppi)
	
	return [monomz,monochg,monointen,monoppi]

def define_charge(mz0_ext,inten0_ext,pmz,maxcharge,unitdiff,ptol):
	
	charge=0
	left=pmz-unitdiff-ptol*pmz*1e-6
	right=pmz+unitdiff+ptol*pmz*1e-6
	pos=nonzero((mz0_ext>=left) & (mz0_ext<=right))[0]
	mz0_n=mz0_ext[pos]
	inten0_n=inten0_ext[pos]
	
	ix = argsort(inten0_n)
	r1 = argsort(ix)
	
	# get charge via m/z on the right with the highest intensity
	for rno in range(len(mz0_n)-1,-1,-1):
		ino=nonzero(r1==rno)[0][0]
		cmz=mz0_n[ino]
		# (1) on the left of pmz, continue
		if cmz<=pmz:
			continue
		# (2)
		diff=1/abs(cmz-pmz)
		# (3)
		round_diff=round(diff)
		# (4)
		if round_diff==0:
			continue
		# (5)
		if round_diff>maxcharge:
			continue
		# (6)
		var=abs(abs(cmz-pmz)-(unitdiff/round_diff))
		# (7)
		if var>ptol*pmz*1e-6:
			continue
		# (8)
		charge=round_diff
		break
	
	# if charge=0, get charge via m/z on the left with the highest intensity
	if charge==0:
		for rno in range(len(mz0_n)-1,-1,-1):
			ino=nonzero(r1==rno)[0][0]
			cmz=mz0_n[ino]
			# (1) on the right of pmz, continue
			if cmz>=pmz:
				continue
			# (2)
			diff=1/abs(cmz-pmz)
			# (3)
			round_diff=round(diff)
			# (4)
			if round_diff==0:
				continue
			# (5)
			if round_diff>maxcharge:
				continue
			# (6)
			var=abs(abs(cmz-pmz)-(unitdiff/round_diff))
			# (7)
			if var>ptol*pmz*1e-6:
				continue
			# (8)
			charge=round_diff
			break
	
	# if charge=1, get charge via mz on the left with the highest intensity
	if charge==1:
		for rno in range(len(mz0_n)-1,-1,-1):
			ino=nonzero(r1==rno)[0][0]
			cmz=mz0_n[ino]
			# (1) on the right of pmz, continue
			if cmz>=pmz:
				continue
			# (2)
			diff=1/abs(cmz-pmz)
			# (3)
			round_diff=round(diff)
			# (4)
			if round_diff==0:
				continue
			# (5)
			if round_diff>maxcharge:
				continue
			# (6)
			var=abs(abs(cmz-pmz)-(unitdiff/round_diff))
			# (7)
			if var>ptol*pmz*1e-6:
				continue
			# (8)
			charge=round_diff
			# (9)
			if charge==1:
				continue
			break
	
	return charge

def get_isotopic_distribution(mz):
	
	loop = 0
	if mz<1500:
		loop = 1
	elif mz<3000:
		loop = 2
	elif mz<4500:
		loop = 4
	elif mz<6000:
		loop = 6
	
	return loop

def get_intensity_ratio(mz,loop):
	
	ratio=2
	if mz<1500:# if the mass <1500, there is no isotopic peaks preceeding the monoisotope
		ratio = 0.3
	elif mz<3000 and loop<=2:
		if loop == 1:
			ratio = 0.2
		elif loop <= 2:
			ratio = 0.1
	elif mz<4500 and loop<=4:
		if loop == 1:
			ratio = 0.2
		elif loop <= 4:
			ratio = 0.1
	elif mz<6000 and loop<=6:
		if loop == 1:
			ratio = 0.2
		elif loop <= 6:
			ratio = 0.1
	
	return ratio

def deisotope(mz0_ext,inten0_ext,mz0,inten0,pno,charge,unitdiff,ptol,procflag):
	
	pmz=mz0[pno]
	pinten=inten0[pno]
	
	mono_mz=pmz
	mono_mz_array=[]
	
	search_loop=1
	flag=1
	mz_hash_saved=pinten
	
	# save cluster peaks +1, +2, etc
	while search_loop>=1 and flag==1:
		flag=0
		left=pmz+unitdiff/charge*search_loop-ptol*pmz*1e-6
		right=pmz+unitdiff/charge*search_loop+ptol*pmz*1e-6
		#pos=nonzero((mz0>=left) & (mz0<=right) & (procflag==0))[0]# no filtering
		pos0=nonzero((mz0>=left) & (mz0<=right))[0]
		ix0=nonzero(procflag[pos0]==0)[0]
		pos=pos0[ix0]
		if len(pos)>0:
			flag=1
			idx=argmax(inten0[pos])
			mz_hash_saved+=inten0[pos[idx]]
			search_loop+=1
			procflag[pno]=1
			procflag[pos[idx]]=1
	
	# save cluster peaks -1
	left=pmz-unitdiff/charge-ptol*pmz*1e-6
	right=pmz-unitdiff/charge+ptol*pmz*1e-6
	#pos=nonzero((mz0>=left) & (mz0<=right) & (procflag==0))[0]# no filtering
	pos0=nonzero((mz0>=left) & (mz0<=right))[0]
	ix0=nonzero(procflag[pos0]==0)[0]
	pos=pos0[ix0]
	if len(pos)>0:
		idx=argmax(inten0[pos])
		mz_hash_saved+=inten0[pos[idx]]
		procflag[pno]=1
		procflag[pos[idx]]=1
	
	# mono_mz
	left=pmz-unitdiff-ptol*pmz*1e-6
	right=pmz+unitdiff+ptol*pmz*1e-6
	pos=nonzero((mz0_ext>=left) & (mz0_ext<=right))[0]
	mz0_n=mz0_ext[pos]
	inten0_n=inten0_ext[pos]
	
	M_zH=pmz*charge
	max_loop=get_isotopic_distribution(M_zH)
	for search_loop in range(0,max_loop+1):
		intensity_ratio=get_intensity_ratio(M_zH,search_loop)
		left=pmz-unitdiff/charge*search_loop-ptol*pmz*1e-6
		right=pmz-unitdiff/charge*search_loop+ptol*pmz*1e-6
		pos=nonzero((mz0_n>=left) & (mz0_n<=right))[0]
		if len(pos)>0:
			idx_set=[]
			for idx in range(0,len(pos)):
				if inten0_n[pos[idx]]/pinten>intensity_ratio:
					idx_set.append(idx)
			if len(idx_set)==1:
				idx=idx_set[0]
				mono_mz=mz0_n[pos[idx]]
				mono_mz_array.append(mono_mz)
			elif len(idx_set)>1:
				ix=argmax(inten0_n[pos[idx_set]])
				mono_mz=mz0_n[pos[idx_set[ix]]]
				mono_mz_array.append(mono_mz)
	
	return [mono_mz,mono_mz_array,mz_hash_saved,procflag]

def changeMH_folder(monomz,monochg,monointen,monoppi,params):
	
	max_num_ppi = params['max_num_ppi']# 0 = disable; 1-10 = max precursor ions selected for mixed MS2 search
	percentage_ppi = params['percentage_ppi']# minimal percentage of precursor peak intensity (ppi) when max_num_ppi = 0
	
	if max_num_ppi<0:
		ix = argsort(monoppi)[::-1]
		r1 = argsort(ix)
		pos=nonzero(array(r1)<=abs(max_num_ppi)-1)[0]
	elif max_num_ppi>0:
		ix = argsort(monoppi)[::-1]
		r1 = argsort(ix)
		pos=nonzero(array(r1)<=max_num_ppi-1)[0]
	else:#max_num_ppi==0:
		pos=nonzero(array(monoppi)*100>percentage_ppi)[0]
	
	if 0==len(pos):
		monomz = []
		monochg = []
		monointen = []
		monoppi = []
	else:
		monomz = list(array(monomz)[pos])
		monochg = list(array(monochg)[pos])
		monointen = list(array(monointen)[pos])
		monoppi = list(array(monoppi)[pos])
	
	return [monomz,monochg,monointen,monoppi]

def GetClusters4Chg1(MS1_index,MS1_peaks,index,params,all_tIPV,ms1_mode,nrunning_mode,cur_premz,cur_ms1pos):
	
	ptola=params['intrascanppm']
	ptole=params['interscanppm']
	unitdiff=params['unitdiff']
	charge_list = params['charge_considered'] # a list of all charges to be considered
	ion_shift = params['monoisotopic_ion_shift'] # a list of monoisotopic ion mass shift to be considered
	left_ion_shift = abs(min(ion_shift))
	max_num_ppi = params['max_num_ppi'] # 0 = disable; 1-10 = max precursor ions selected for mixed MS2 search
	percentage_ppi = params['percentage_ppi'] # minimal percentage of precursor peak intensity (ppi) when max_num_ppi = 0
	simple_process = params['simple_process'] # if 'deisotoping_method=1 or 2', 1: simple process with 2 steps; 2: complex process with 6 steps; if 'deisotoping_method=0', 1: charge +2; 2: charges +2, +3
	chg_range1=arange(1,1+1)#only 1
	
	monomz = []
	monochg = []
	monointen = []
	monoppi = []
	monogrp = []
	
	if simple_process!=1:
		#complex
		cmz=[]
		grp=[]
		[sel_pnum,sel_pos,sel_mz,sel_inten] = SelectPeaks(MS1_index,MS1_peaks,index,params,ms1_mode,cur_premz,cur_ms1pos)
		[sel_pnum,sel_pos,sel_mz,sel_inten] = CheckAdjMS1scans(params,ms1_mode,cur_premz,cur_ms1pos,sel_pnum,sel_pos,sel_mz,sel_inten,1,'')
		[mz,inten] = GetRefMS1(sel_pnum,sel_mz,sel_inten,ptole)
		[isonum,isomz,isointen,chg] = GetAllClustersOnSingleMS1(mz,inten,charge_list,unitdiff,ptola*2,1)
		inum=size(isomz,0)
		if 0==inum:
			return [monomz,monochg,monointen,monoppi]
		
		[isonum,isomz,isointen,chg]=FilterByPPI(isonum,isomz,isointen,chg,0,cur_premz,ptola*2)
		if size(isomz,0)>0:
			[cmz,grp]=Get_cmz_from_clusters_curpremz(isonum,isomz,isointen,chg,left_ion_shift,cur_premz,ptola*2)
		else:
			cmz.append(cur_premz)
			grp.append(100)
	else:
		# #simple
		cmz=[]
		grp=[]
		[sel_pnum,sel_pos,sel_mz,sel_inten] = SelectPeaks(MS1_index,MS1_peaks,index,params,ms1_mode,cur_premz,cur_ms1pos)
		[mz,inten] = GetOneMS1(sel_pnum,sel_pos,sel_mz,sel_inten,ms1_mode,cur_ms1pos,0)
		cmz.append(cur_premz)
		grp.append(100)
	
	gno=0
	for cno in range(0,len(cmz)):
		for dno in range(0,len(chg_range1)):
			gno=grp[cno]+dno+1
			for eno in range(0,len(ion_shift)):
				monomz.append(cmz[cno]+ion_shift[eno]*unitdiff/chg_range1[dno])
				monochg.append(chg_range1[dno])
				monointen.append(100)
				monoppi.append(0.01)
				monogrp.append(gno)
	
	[isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp] = DeisotopeMS1Features(monomz,monochg,monogrp,params,mz,inten,ptola*2,all_tIPV)
	inum = size(isomz,0)
	if 0==inum:
		return [monomz,monochg,monointen,monoppi]
	[isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp] = FilterByISO(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,ptola*2)
	[monomz,monochg,monointen,monoppi] = GetTopPPI(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,max_num_ppi,percentage_ppi)
	return [monomz,monochg,monointen,monoppi]

def GetClusters(MS1_index,MS1_peaks,index,params,all_tIPV,ms1_mode,nrunning_mode,cur_premz,cur_ms1pos):
	# Get all isotopic cluster candidates in the window
	
	ptola=params['intrascanppm']
	ptole=params['interscanppm']
	unitdiff=params['unitdiff']
	nHalfWidthL=params['half_width_left']
	nHalfWidthR=params['half_width_right']
	topN = params['precursor_ion_considered'] # N (1-4): go back to MS1, start from the strongest peak of top-N clusters in the isolation window
	charge_list = params['charge_considered'] # a list of all charges to be considered
	ion_shift = params['monoisotopic_ion_shift'] # a list of monoisotopic ion mass shift to be considered
	left_ion_shift = abs(min(ion_shift))
	maxcharge = max(charge_list)
	max_num_ppi = params['max_num_ppi'] # 0 = disable; 1-10 = max precursor ions selected for mixed MS2 search
	percentage_ppi = params['percentage_ppi'] # minimal percentage of precursor peak intensity (ppi) when max_num_ppi = 0
	simple_process = params['simple_process'] # if 'deisotoping_method=1 or 2', 1: simple process with 2 steps; 2: complex process with 6 steps; if 'deisotoping_method=0', 1: charge +2; 2: charges +2, +3
	chg_range2=array([i for i in charge_list if i!=1])#start from 2
	
	monomz = []
	monochg = []
	monointen = []
	monoppi = []
	monogrp = []
	
	# via_running_mode
	# 1: brute force; 2: pattern matching via entire envelope; 3: step-wise via adjacent peaks; 4: debug
	if nrunning_mode==1 or nrunning_mode==2:
		if simple_process!=1:
			#complex
			cmz=[]
			grp=[]
			[sel_pnum,sel_pos,sel_mz,sel_inten] = SelectPeaks(MS1_index,MS1_peaks,index,params,ms1_mode,cur_premz,cur_ms1pos)
			[sel_pnum,sel_pos,sel_mz,sel_inten] = CheckAdjMS1scans(params,ms1_mode,cur_premz,cur_ms1pos,sel_pnum,sel_pos,sel_mz,sel_inten,2,'')
			[mz,inten] = GetRefMS1(sel_pnum,sel_mz,sel_inten,ptole)
			[isonum,isomz,isointen,chg] = GetAllClustersOnSingleMS1(mz,inten,charge_list,unitdiff,ptola*2,2)
			inum=size(isomz,0)
			if 0==inum:
				# ########if no >=+2 charges, check charge +1 only########
				[monomz,monochg,monointen,monoppi]=GetClusters4Chg1(MS1_index,MS1_peaks,index,params,all_tIPV,ms1_mode,nrunning_mode,cur_premz,cur_ms1pos)
				return [monomz,monochg,monointen,monoppi]
			
			[isonum,isomz,isointen,chg]=FilterByPPI(isonum,isomz,isointen,chg,topN,cur_premz,ptola*2)
			if size(isomz,0)>0:
				[cmz,grp]=Get_cmz_from_clusters_curpremz(isonum,isomz,isointen,chg,left_ion_shift,cur_premz,ptola*2)
			else:
				cmz.append(cur_premz)
				grp.append(100)
		else:
			# #simple
			cmz=[]
			grp=[]
			[sel_pnum,sel_pos,sel_mz,sel_inten] = SelectPeaks(MS1_index,MS1_peaks,index,params,ms1_mode,cur_premz,cur_ms1pos)
			[mz,inten] = GetOneMS1(sel_pnum,sel_pos,sel_mz,sel_inten,ms1_mode,cur_ms1pos,0)
			cmz.append(cur_premz)
			grp.append(100)
		
		gno=0
		for cno in range(0,len(cmz)):
			for dno in range(0,len(chg_range2)):
				gno=grp[cno]+dno+1
				for eno in range(0,len(ion_shift)):
					monomz.append(cmz[cno]+ion_shift[eno]*unitdiff/chg_range2[dno])
					monochg.append(chg_range2[dno])
					monointen.append(100)
					monoppi.append(0.01)
					monogrp.append(gno)
		
		if nrunning_mode==1:
			return [monomz,monochg,monointen,monoppi]
		elif nrunning_mode==2:
			[isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp] = DeisotopeMS1Features(monomz,monochg,monogrp,params,mz,inten,ptola*2,all_tIPV)
			inum = size(isomz,0)
			if 0==inum:
				# ########if no >=+2 charges, check charge +1 only########
				[monomz,monochg,monointen,monoppi]=GetClusters4Chg1(MS1_index,MS1_peaks,index,params,all_tIPV,ms1_mode,nrunning_mode,cur_premz,cur_ms1pos)
				return [monomz,monochg,monointen,monoppi]
			[isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp] = FilterByISO(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,ptola*2)
			[monomz,monochg,monointen,monoppi] = GetTopPPI(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,max_num_ppi,percentage_ppi)
			return [monomz,monochg,monointen,monoppi]
	elif nrunning_mode==3:
		if params['ncorrect']==0:# if uncorrected, double params['intrascanppm']
			ptola = ptola*2
		[sel_pnum,sel_pos,sel_mz,sel_inten] = SelectPeaks(MS1_index,MS1_peaks,index,params,ms1_mode,cur_premz,cur_ms1pos)
		[mz0_ext,inten0_ext] = GetOneMS1(sel_pnum,sel_pos,sel_mz,sel_inten,ms1_mode,cur_ms1pos,0)
		
		left_lim=cur_premz-nHalfWidthL
		right_lim=cur_premz+nHalfWidthR
		pos=nonzero((mz0_ext>=left_lim) & (mz0_ext<=right_lim))[0]
		if len(pos)==0:
			return [monomz,monochg,monointen,monoppi]
		mz0=mz0_ext[pos]# from extend to isolation window
		inten0=inten0_ext[pos]
		
		# r1: rank of inten0 by ascending order (the last is the highest)
		ix = argsort(inten0)
		r1 = argsort(ix)
		
		# decharge through the rank (the highest the first)
		procflag=zeros((len(mz0),1))
		for rno in range(len(mz0)-1,-1,-1):
			ino=nonzero(r1==rno)[0][0]
			if 1==procflag[ino]:
				continue
			
			charge=define_charge(mz0_ext,inten0_ext,mz0[ino],maxcharge,unitdiff,ptola)
			
			if charge!=0:
				[mono_mz,mono_mz_array,mz_hash_saved,procflag]=deisotope(mz0_ext,inten0_ext,mz0,inten0,ino,charge,unitdiff,ptola,procflag)
				monomz.append(mono_mz)
				monochg.append(charge)
				monointen.append(mz_hash_saved)
			else:
				monomz.append(mz0[ino])
				monochg.append(0)
				monointen.append(inten0[ino])
				procflag[ino]=1
		
		# PPI
		if len(monointen)==0:
			return [monomz,monochg,monointen,monoppi]
		intensity_sum = sum(monointen)
		for ino in range(0,len(monointen)):
			monoppi.append(monointen[ino]/intensity_sum)
		
		# changeMH folder
		[monomz,monochg,monointen,monoppi]=changeMH_folder(monomz,monochg,monointen,monoppi,params)
		
		return [monomz,monochg,monointen,monoppi]

def Get_cmz_from_clusters_curpremz(isonum,isomz,isointen,chg,left_ion_shift,cur_premz,ptol):
	# Get cmz from clusters cur_premz and get rid of the redundant
	
	# check if cur_premz in clusters
	inum=size(isomz,0)
	c_ptol=ptol*cur_premz*1e-6
	ninclude_cur_premz=0
	flag_ppm=[]# how many mz in the same cluster
	idx_ref=[]# where is mz from (the cur_premz or the highest peak)
	for ino in range(inum):
		c_isomz=isomz[ino,0:isonum[ino]]
		IX=nonzero((cur_premz-c_ptol<=c_isomz) & (c_isomz<=cur_premz+c_ptol) )[0]
		if len(IX)>0:
			ninclude_cur_premz=1
			isomz[ino,IX[0]]=cur_premz
			idx_ref.append(IX[0])
			if IX[0]<=left_ion_shift:# in the range of left ion shift
				flag_ppm.append(1)
			else:# out the range of left ion shift, use two m/z
				flag_ppm.append(2)
		else:
			c_isointen=isointen[ino,0:isonum[ino]]
			idx=argmax(c_isointen)
			idx_ref.append(idx)
			if idx<=left_ion_shift:# in the range of left ion shift
				flag_ppm.append(1)
			else:# out the range of left ion shift, use two m/z
				flag_ppm.append(2)
	flag_ppm=array(flag_ppm)
	idx_ref=array(idx_ref)
	
	x=flag_ppm.argsort()[::-1]# sort flag_ppm in reverse order
	isonum=isonum[x]
	isomz=isomz[x,:]
	isointen=isointen[x,:]
	chg=chg[x]
	flag_ppm=flag_ppm[x]
	idx_ref=idx_ref[x]
	
	cmz=[]
	grp=[]
	gno=0
	for ino in range(inum):
		c_isomz=isomz[ino,0:isonum[ino]]
		if flag_ppm[ino]==1:# in the range of left ion shift
			cmz.append( c_isomz[idx_ref[ino]] )
			gno=gno+1
			grp.append(gno*100)
		else:# out the range of left ion shift, use two m/z
			cmz.append( c_isomz[idx_ref[ino]] )
			cmz.append( c_isomz[0] )
			gno=gno+1
			grp.append(gno*100)
			grp.append(gno*100)
	# cur_premz not in clusters, add cur_premz
	if ninclude_cur_premz==0:
		cmz.append(cur_premz)
		gno=gno+1
		grp.append(gno*100)
	
	# get rid of the redundant
	cmz_array=array(cmz)
	grp_array=array(grp)
	flag=array([1]*len(cmz_array))
	if len(cmz_array)>1:
		for ino in range(0,len(cmz_array)-1):
			mz_i=cmz_array[ino]
			c_ptol = ptol*mz_i*1e-6
			for jno in range(ino+1,len(cmz_array)):
				mz_j=cmz_array[jno]
				if abs(mz_j-mz_i)<=c_ptol:
					flag[jno]=0
		ix=nonzero(flag==1)[0]
		cmz_array=cmz_array[ix]
		grp_array=grp_array[ix]
	cmz=cmz_array.tolist()
	grp=grp_array.tolist()
	
	return [cmz,grp]

def GetClusters_psm(MS1_index,MS1_peaks,index,params,all_tIPV,ms1_mode,nrunning_mode,cur_premz,cur_ms1pos,cur_ms2scan,logfile):
	# Get all isotopic cluster candidates in the window
	
	ptola=params['intrascanppm']
	ptole=params['interscanppm']
	unitdiff=params['unitdiff']
	nHalfWidthL=params['half_width_left']
	nHalfWidthR=params['half_width_right']
	topN = params['precursor_ion_considered'] # N (1-4): go back to MS1, start from the strongest peak of top-N clusters in the isolation window
	charge_list = params['charge_considered'] # a list of all charges to be considered
	ion_shift = params['monoisotopic_ion_shift'] # a list of monoisotopic ion mass shift to be considered
	left_ion_shift = abs(min(ion_shift))
	max_num_ppi = params['max_num_ppi'] # 0 = disable; 1-10 = max precursor ions selected for mixed MS2 search
	percentage_ppi = params['percentage_ppi'] # minimal percentage of precursor peak intensity (ppi) when max_num_ppi = 0
	chg_range2=array([i for i in charge_list if i!=1])#start from 2
	psm_scan=params['target_ms2scan']
	psm=params['target_psm']
	
	monomz = []
	monochg = []
	monointen = []
	monoppi = []
	monogrp = []
	
	# psm
	if nrunning_mode==4:
		if len(psm)>0:
			psm_scan=array(psm_scan)
			X=nonzero( psm_scan==cur_ms2scan )[0]
			if len(X)==0:
				return [monomz,monochg,monointen,monoppi]
			
			for t in range(0,len(X)):
				i=X[t]
				ntype=check_pre_type_simple(psm[i][2],psm[i][1],cur_premz)
				logfile.write('%d\t%d\t%f\t%s\t%s\t%f\t%s\t%d\n' % (psm[i][0],psm[i][1],psm[i][2],psm[i][3],psm[i][4],psm[i][5],psm[i][6],ntype))
			
			return [monomz,monochg,monointen,monoppi]
		else:
			cmz=[]
			grp=[]
			[sel_pnum,sel_pos,sel_mz,sel_inten] = SelectPeaks(MS1_index,MS1_peaks,index,params,ms1_mode,cur_premz,cur_ms1pos)
			[sel_pnum,sel_pos,sel_mz,sel_inten] = CheckAdjMS1scans(params,ms1_mode,cur_premz,cur_ms1pos,sel_pnum,sel_pos,sel_mz,sel_inten,2,logfile)
			[mz,inten] = GetRefMS1(sel_pnum,sel_mz,sel_inten,ptole)
			[isonum,isomz,isointen,chg] = GetAllClustersOnSingleMS1(mz,inten,charge_list,unitdiff,ptola*2,2)
			
			[isonum,isomz,isointen,chg]=FilterByPPI(isonum,isomz,isointen,chg,topN,cur_premz,ptola*2)
			if size(isomz,0)>0:
				[cmz,grp]=Get_cmz_from_clusters_curpremz(isonum,isomz,isointen,chg,left_ion_shift,cur_premz,ptola*2)
			else:
				cmz.append(cur_premz)
				grp.append(100)
			
			gno=0
			for cno in range(0,len(cmz)):
				for dno in range(0,len(chg_range2)):
					gno=grp[cno]+dno+1
					for eno in range(0,len(ion_shift)):
						monomz.append(cmz[cno]+ion_shift[eno]*unitdiff/chg_range2[dno])
						monochg.append(chg_range2[dno])
						monointen.append(100)
						monoppi.append(0.01)
						monogrp.append(gno)
			
			nchoice=2
			
			if nchoice==1:
				# 1
				logfile.write('----------------GetClusters----------------\n')
				logfile.write('npre	premz	prech	pretype	prePPI	preDscore	predelta_Dscore\n')
				logfile.write('%d' % (len(cmz)))
				for i in range(0,len(cmz)):
					logfile.write('\t%f\t%d\t%d\t%f\t%f\t%f\n' % (cmz[i],0,0,0,0,0))
				logfile.write('\n')
				
				return [monomz,monochg,monointen,monoppi]
			else:
				# 2
				[isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp] = DeisotopeMS1Features(monomz,monochg,monogrp,params,mz,inten,ptola*2,all_tIPV)
				inum = size(isomz,0)
				if 0==inum:
					monomz = []
					monochg = []
					monointen = []
					monoppi = []
					return [monomz,monochg,monointen,monoppi]
				[isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp] = FilterByISO(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,ptola*2)
				[monomz,monochg,monointen,monoppi] = GetTopPPI(isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp,max_num_ppi,percentage_ppi)
				
				inum = size(isomz,0)
				logfile.write('----------------GetClusters----------------\n')
				logfile.write('npre	premz	prech	pretype	prePPI	preDscore	predelta_Dscore\n')
				logfile.write('%d' % (inum))
				for i in range(0,inum):
					logfile.write('\t%f\t%d\t%d\t%f\t%f\t%f\n' % (isomz[i,0],chg[i],0,suminten[i]/sum(suminten),Dscore[i],delta_Dscore[i]))
				logfile.write('\n')
				logfile.write('\n')
				logfile.write('----------------GetTopPPI----------------\n')
				logfile.write('	monomz	monochg	monointen	monoppi\n')
				for i in range(0,len(monomz)):
					logfile.write('\t%f\t%d\t%f\t%f\n' % (monomz[i],monochg[i],monointen[i],monoppi[i]))
				logfile.write('\n')
				
				return [monomz,monochg,monointen,monoppi]
	
	return [monomz,monochg,monointen,monoppi]

def check_pre_type_simple(psm_mz,psm_ch,cur_premz):
	
	ntype=9
	
	C=1.00335
	ptol=20
	L=[0,-1,-2,-3,-4,-5,-6,1,2,3,4,5,6]
	c_ptol=ptol*psm_mz*1e-6
	for k in range(0,len(L)):
		c_mz=cur_premz
		c_ch=psm_ch
		check_list=c_mz+L[k]*C/c_ch
		check_list=array(check_list)
		IX=nonzero((psm_mz-c_ptol<=check_list) & (check_list<=psm_mz+c_ptol) )[0]
		if len(IX)>0:
			ntype=L[k]
			break
	
	return ntype

def GetMonoPeaks(c_mzXML_fullfile,jump_params,cur_nbatch,total_nbatch,ncorrect):
	# get monoisotopics
	
	# load ms1 and ms2
	[tpath,tname] = os.path.split(c_mzXML_fullfile)
	[fname,ext] = os.path.splitext(tname)
	t_mzXML_fullfile = os.path.join(tpath,fname,tname)
	
	MS2_monofile = Change_ext(t_mzXML_fullfile,'_Mono%d.npz' % (cur_nbatch))
	if os.path.isfile(MS2_monofile):
		return
	
	MS1_scanfile = Change_ext(t_mzXML_fullfile,'_MS1scans.npz')
	MS1_peakfile = Change_ext(t_mzXML_fullfile,'_MS1peaks.npz')
	MS2_scanfile = Change_ext(t_mzXML_fullfile,'_MS2scans.npz')
	MS2_peakfile = Change_ext(t_mzXML_fullfile,'_MS2peaks.npz')
	
	M1scan=load(MS1_scanfile)
	M1peak=load(MS1_peakfile)
	MS1_index=M1scan['MS1_index']
	MS1_peaks=M1peak['MS1_peaks']
	num_MS1=MS1_index.shape[0]
	M2scan=load(MS2_scanfile)
	MS2_index=M2scan['MS2_index']
	num_MS2=MS2_index.shape[0]
	#print([num_MS1,num_MS2])
	index0=[1]+list(MS1_index[0:num_MS1,2])
	index=[int(i) for i in index0]
	#print([index[0],index[1],len(index)])
	
	# get cur_num_MS2
	cur_MS2POS=get_cur_MS2POS(num_MS2,cur_nbatch,total_nbatch)
	cur_num_MS2=len(cur_MS2POS)
	
	maxprenum=1000
	ms2scan = zeros((cur_num_MS2,1))
	mononum = zeros((cur_num_MS2,1))
	monomz = zeros((cur_num_MS2,maxprenum))
	monochg = zeros((cur_num_MS2,maxprenum))
	monointen = zeros((cur_num_MS2,maxprenum))
	monoppi = zeros((cur_num_MS2,maxprenum))
	
	# get internal parameters
	params = Get_params(jump_params)
	params['ncorrect'] = ncorrect # ###add ncorrect to params###
	nTMT = len(params['tmt'])	# TMT reporter ion mass table
	all_tIPV=Getall_tIPV(nTMT)
	ms1_mode = Get_ms1_mode(MS1_peaks,index)
	nrunning_mode=params['deisotoping_method']
	ntarget_ms2scan=params['target_ms2scan']
	psm=params['target_psm']
	simple_process = params['simple_process'] # if 'deisotoping_method=1 or 2', 1: simple process with 2 steps; 2: complex process with 6 steps; if 'deisotoping_method=0', 1: charge +2; 2: charges +2, +3
	
	if nrunning_mode==0:
		# mass correction only
		ms2scan = zeros((cur_num_MS2,1))
		mononum = zeros((cur_num_MS2,1))
		monomz = zeros((cur_num_MS2,1))
		monochg = zeros((cur_num_MS2,1))
		monointen = zeros((cur_num_MS2,1))
		monoppi = zeros((cur_num_MS2,1))
		for cno in range(0,cur_num_MS2):
			# Each MS2 scan
			sno=cur_MS2POS[cno]
			
			# MS1_index=zeros((MaxMSNum,5)) # MS1 scan, MS1 rt, MS1 peak num, baseline, massshift
			# MS2_index=zeros((MaxMSNum,8)) # MS1 scan,MS1 rt,MS2 scan,m/z,z,Fragtype,MS2 peak num,massshift
			cur_ms2scan=int(MS2_index[sno,2])
			cur_premz=MS2_index[sno,3]
			cur_prech=int(MS2_index[sno,4])
			
			ms2scan[cno]=cur_ms2scan
			mononum[cno]=1
			monomz[cno]=cur_premz
			monochg[cno]=cur_prech
			monointen[cno]=1e6
			monoppi[cno]=1
			if simple_process!=1:
				monochg[cno]=0
		
		savez(MS2_monofile,ms2scan=ms2scan,mononum=mononum,monomz=monomz,monochg=monochg,monointen=monointen,monoppi=monoppi)
		return
	
	if nrunning_mode==4:
		# MS2_logfile
		MS2_logfile = Change_ext(t_mzXML_fullfile,'_log%d.xls' % (cur_nbatch))
		logfile = open(MS2_logfile,'w')
		if cur_nbatch==1 and len(psm)>0:
			logfile.write('scan	ch	mz	sq	mod	score	ac	type\n')
	
	# Flowchart
	#print('ms2 scans(%d):' % (cur_num_MS2))
	for cno in range(0,cur_num_MS2):
		# Each MS2 scan
		#print(cno+1, end="\r")
		sno=cur_MS2POS[cno]
		
		# MS1_index=zeros((MaxMSNum,5)) # MS1 scan, MS1 rt, MS1 peak num, baseline, massshift
		# MS2_index=zeros((MaxMSNum,8)) # MS1 scan,MS1 rt,MS2 scan,m/z,z,Fragtype,MS2 peak num,massshift
		cur_ms2scan=int(MS2_index[sno,2])
		cur_premz=MS2_index[sno,3]
		cur_prech=int(MS2_index[sno,4])
		cur_ms1scan=int(MS2_index[sno,0])
		cur_ms1pos=nonzero(MS1_index[:,0]==cur_ms1scan)[0][0]
		#print([cur_ms2scan,cur_premz,cur_prech,cur_ms1scan,cur_ms1pos])
		
		if 0==cur_ms1pos.size:
			c_monomz = []
			c_monochg = []
			c_monointen = []
			c_monoppi = []
		else:
			if nrunning_mode!=4:
				[c_monomz,c_monochg,c_monointen,c_monoppi]=GetClusters(MS1_index,MS1_peaks,index,params,all_tIPV,ms1_mode,nrunning_mode,cur_premz,cur_ms1pos)
			else:#nrunning_mode=4, debug mode
				if cur_ms2scan in ntarget_ms2scan:
					[c_monomz,c_monochg,c_monointen,c_monoppi]=GetClusters_psm(MS1_index,MS1_peaks,index,params,all_tIPV,ms1_mode,nrunning_mode,cur_premz,cur_ms1pos,cur_ms2scan,logfile)
					if len(psm)==0:
						logfile.write('----------------GetMonoPeaks----------------\n')
						logfile.write('cur_premz:\n%.2f-%.2f\ncur_ms2rt:\n%.2f\ncur_ms2scan:\n%d\n' % (cur_premz-2,cur_premz+2,MS2_index[sno,1],cur_ms2scan))
				else:
					c_monomz = []
					c_monochg = []
					c_monointen = []
					c_monoppi = []
				
		ms2scan[cno]=cur_ms2scan
		inum = len(c_monomz)
		if inum>0:
			mononum[cno]=inum
			monomz[cno,0:inum]=c_monomz
			monochg[cno,0:inum]=c_monochg
			monointen[cno,0:inum]=c_monointen
			monoppi[cno,0:inum]=c_monoppi
	#print('')
	
	if nrunning_mode==4:
		# MS2_logfile
		logfile.close()
	
	mno=int(max(mononum))
	if mno<maxprenum:
		monomz=monomz[:,0:mno]
		monochg=monochg[:,0:mno]
		monointen=monointen[:,0:mno]
		monoppi=monoppi[:,0:mno]
	if nrunning_mode!=4:
		savez(MS2_monofile,ms2scan=ms2scan,mononum=mononum,monomz=monomz,monochg=monochg,monointen=monointen,monoppi=monoppi)
	
	return
