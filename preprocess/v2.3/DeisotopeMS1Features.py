from numpy import *

def DeisotopeMS1Features(monomz,monochg,monogrp,params,mz,inten,ptol,all_tIPV):
	# DeisotopeMS1Features
	
	# params
	tags_list = params['number_TMT_tags'] # a list of TMT tags to be considered
	IsoIntenCutoff = params['isotopic_pattern_cutoff'] # percentage threshold (relative to the strongest ion) that a theoretical ion to be considered
	delta_DscoreCutoff = params['delta_Dscore'] # threshold within which multiple precursors may be reported; delta_Dscore = (Dscore_max - Dscore_current) / Dscore_max
	unitdiff=params['unitdiff']
	Tsim=params['tsim']
	
	mz_len=len(monomz)
	
	# isotopic clusters
	isonum=array(mz_len*[0])
	isomz=zeros((mz_len,60))
	isointen=zeros((mz_len,60))
	chg=array(mz_len*[0])
	
	# features
	Dscore=array([0.0]*mz_len)
	delta_Dscore=array([0.0]*mz_len)
	suminten=array([0.0]*mz_len)
	isogrp=array(monogrp)
	flag=array([0]*mz_len)
	
	#Dscore
	for ino in range(mz_len):
		c_pmz=monomz[ino]
		c_chg=monochg[ino]
		
		# experiment clusters
		c_mz = [c_pmz-unitdiff/c_chg,c_pmz,c_pmz+unitdiff/c_chg,c_pmz+2*unitdiff/c_chg,c_pmz+3*unitdiff/c_chg,c_pmz+4*unitdiff/c_chg,c_pmz+5*unitdiff/c_chg]
		[e_mz,e_inten] = AlignMS1(mz,inten,c_mz,ptol)
		
		# features (Dscore, suminten)
		[Dscore[ino],tm]=get_sim(c_pmz,c_chg,e_inten,tags_list,IsoIntenCutoff,all_tIPV)
		e_mz=e_mz[1:tm]
		e_inten=e_inten[1:tm]
		suminten[ino]=sum(e_inten)
		
		# isotopic clusters
		isonum[ino]=len(e_mz)
		isomz[ino,0:len(e_mz)]=e_mz
		isointen[ino,0:len(e_mz)]=e_inten
		chg[ino]=c_chg
	
	#delta_Dscore, flag
	max_Dscore=max(Dscore)+1e-10
	for ino in range(mz_len):
		delta_Dscore[ino]=(max_Dscore-Dscore[ino])/max_Dscore
		if delta_Dscore[ino]<delta_DscoreCutoff and Dscore[ino]>=Tsim:
			flag[ino]=1
	
	#flag
	idx=nonzero(flag==1)[0]
	if len(idx)>0:
		isonum = isonum[idx]
		isomz = isomz[idx,:]
		isointen = isointen[idx,:]
		chg = chg[idx]
		
		Dscore=Dscore[idx]
		delta_Dscore=delta_Dscore[idx]
		suminten=suminten[idx]
		isogrp=isogrp[idx]
		suminten=recalculate_suminten(suminten,isogrp)
	else:
		mz_len=0
		isonum=array(mz_len*[0])
		isomz=zeros((mz_len,60))
		isointen=zeros((mz_len,60))
		chg=array(mz_len*[0])
		
		Dscore = array([])
		delta_Dscore = array([])
		suminten = array([])
		isogrp = array([])
	
	return [isonum,isomz,isointen,chg,Dscore,delta_Dscore,suminten,isogrp]

def AlignMS1(mz,inten,mz_n,ptol):
	
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

def get_sim(c_pmz,c_chg,e_inten,tags_list,IsoIntenCutoff,all_tIPV):
	#
	pmass=1.007276
	M = c_pmz*c_chg-c_chg*pmass
	
	if len(all_tIPV)==1:
		TMT_data=0
	else:
		TMT_data=1
	
	if TMT_data==0:
		all_tIPV_0 = all_tIPV["0"]
		tIPV_0=all_tIPV_0[min(len(all_tIPV_0)-1,int(M)-1)]
		
		[Dscore,tm]=get_Dscore(e_inten,tIPV_0,IsoIntenCutoff)
	else:
		scores=[]
		tms=[]
		for ino in range(0,len(tags_list)):
			ntag=tags_list[ino]
			if ntag>4 or ntag<0:
				continue
			
			all_tIPV_0 = all_tIPV[str(ntag)]
			tIPV_0=all_tIPV_0[min(len(all_tIPV_0)-1,int(M)-1)]
			[Dscore0,tm0]=get_Dscore(e_inten,tIPV_0,IsoIntenCutoff)
			scores.append(Dscore0)
			tms.append(tm0)
		
		Dscore=max(scores)
		idx=argmax(scores)
		tm=tms[idx]
	
	return [Dscore,tm]

def get_Dscore(e_inten,tIPV,IsoIntenCutoff):
	
	nlim=len(tIPV)
	t_inten=array([0.0]*nlim)
	for tno in range(0,nlim):
		t_inten[tno] = tIPV[tno]
	
	# cut short clusters by relative intensity
	for tno in range(3,nlim):
		if t_inten[tno]<IsoIntenCutoff:
			tm=tno
			break
	if t_inten[nlim-1]>=IsoIntenCutoff:
		tm=nlim
	e_inten=e_inten[0:tm]
	t_inten=t_inten[0:tm]
	
	Dscore0_global=get_similarity(e_inten,t_inten)
	Dscore0_local=get_similarity(e_inten[0:3],t_inten[0:3])
	Dscore0=0.5*Dscore0_global+0.5*Dscore0_local
	
	return [Dscore0,tm]

def get_similarity(e_intens,t_intens):
	
	if sqrt( sum(e_intens*e_intens)*sum(t_intens*t_intens) )==0 or e_intens[1]==0:
		e_sim=0
	else:
		e_max = e_intens[1]
		e_intens_norm = e_intens/e_max*100
		condi_inten0 = (e_intens_norm[0]<70)#(t_intens[0]==0 and e_intens_norm[0]<10) or (t_intens[0]>0 and e_intens_norm[0]<min([2*t_intens[0],40]))
		condi_inten2 = e_intens_norm[2]>1/3*t_intens[2] and e_intens_norm[2]<3*t_intens[2]
		if condi_inten0 and condi_inten2:
			e_sim = sum(e_intens*t_intens)/sqrt( sum(e_intens*e_intens)*sum(t_intens*t_intens) )
		else:
			e_sim=0
	
	return e_sim

def recalculate_suminten(suminten,isogrp):
	# recalculate suminten in the same group
	grps=unique(isogrp)
	for ino in range(0,len(grps)):
		ix=nonzero(isogrp==grps[ino])[0]
		if len(ix)>1:
			suminten_max_in_grp = max(suminten[ix])+0.0
			suminten_sum_in_grp = sum(suminten[ix])+0.0
			suminten[ix] = suminten[ix]/suminten_sum_in_grp*suminten_max_in_grp
	
	return suminten
