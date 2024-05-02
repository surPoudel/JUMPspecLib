#!/usr/bin/python

import os
import sys
import configparser
from DeisotopeMS1Mono import *

def main():
	# ******start****** #
	
	# params_file
	if len(sys.argv)<2:
		print('please input a parameter file!')
		return
	else:
		params_file = sys.argv[1].strip()
		if False==os.path.isfile(params_file):
			print('please input a correct parameter file!')
			return
	
	# read params_file
	cfg = configparser.ConfigParser()
	cfg.read(params_file)
	
	datapath = cfg['Batch']['datapath']
	c_mzXML_fullfile = datapath+'.mzXML'
	jump_params = cfg['Batch']['jump_params']
	cur_nbatch = int(float(cfg['Batch']['current_batch']))
	total_nbatch = int(float(cfg['Batch']['total_batch']))
	ncorrect = int(float(cfg['Batch']['ncorrect']))
	
	# get monoisotopics on current batch
	GetMonoPeaks(c_mzXML_fullfile,jump_params,cur_nbatch,total_nbatch,ncorrect)
	
	return

if __name__ == '__main__':
	main()
