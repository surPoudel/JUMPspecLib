import pandas as pd 
import re
import os, sys

from os.path import dirname
from TargetDecoy import *
from logFunctions import *

def mergeSearchResults(list_of_searches, filter_col= "JDscore", jdscore_cutoff=0.6):
    all_frames = []
    for row in list_of_searches:
        searchFile = row
        df_s = pd.read_csv(searchFile)
        df_s = df_s.query("{}>{}".format(filter_col, jdscore_cutoff))
        if "deltaRT_postcal" in df_s.columns:
            df_s.rename(columns={df_s.columns[20]:"expRT",df_s.columns[21]:"nPSMs"},inplace=True)
        write_log("Input PSMs for initial filtering based on JDScore for {} is {}\n".format(searchFile,df_s.shape[0]))
        all_frames.append(df_s)
    df_s_filtered_co = pd.concat(all_frames)
    return df_s_filtered_co



