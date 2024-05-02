#!/usr/bin/env python3
import os
import sys
import argparse
from os.path import dirname


def msg(name=None):
    return '''\njump_lib <mode> <arguments>\n    mode is one of "-pp", "-mc", "-d", "-d_merge", "-s", "-pp_s","-f", or "-q"\n\n\n '''


parser = argparse.ArgumentParser(description="JUMP suite to perform spectral libray generation, search, filter and quantification", prog="jump_lib_wrapper",usage=msg())
parser.add_argument('-pp', '--preprocess', action='store_true')
# parser.add_argument('-mc', '--masscal', action='store_true')
parser.add_argument('-d', '--makeLibrary', action='store_true')
parser.add_argument('-d_merge', '--mergeLibraries', action='store_true')
parser.add_argument('-s', '--librarySearch', action='store_true')
parser.add_argument('-pp_s', '--librarySearch_w_pp', action='store_true')
parser.add_argument('-f', '--libraryFilter', action='store_true')
parser.add_argument('-q', '--libraryQuan', action='store_true')

parser.add_argument("parameterfile",help="mode specific parameter file for example for --makeLibrary specLibGen.params")
parser.add_argument("mzXML",help="single or list of mzXML files",nargs='*')



parser.add_argument("--queue",dest="queue", default="queue=standard")
parser.add_argument("--mem",dest="mem", default="mem=20000")

args = parser.parse_args()

# path for the script
source_path_script = sys.argv[0]
#get the path of script
source_path = dirname(source_path_script)
# print (source_path)

paramterFile = args.parameterfile
if args.mzXML:
    inputFile = " ".join(args.mzXML)
    # print (inputFile)


queue = args.queue.split("=")[-1]
mem = args.mem.split("=")[-1]
preprocess = "{}/preprocess/preprocess_lib.py".format(source_path)
libgen = "{}/consensusLibrary/consensusLibraryGeneration.py".format(source_path)
libmerge = "{}/consensusLibrary/DatabaseMerging.py".format(source_path)
libsearch = "{}/JumplibrarySearch/librarySearchMain.py".format(source_path)
libsearch_w_preprocess = "{}/JumplibrarySearch/librarySearchMain_withPreprocess.py".format(source_path)
libfilter = "{}/JumplibraryFilter/Jump_lib_filter_2modes.py".format(source_path)
libquan = "{}/JumplibraryQuan/jump_lib_quan_v0.1.0.py".format(source_path)


# #mode is either -d, or -d_merge or -s or -f or -q

# # print (paramterFile)

# if (len(mode) == 1) and ("d" in mode):
if args.makeLibrary:
    cmd = "python {} {}".format(libgen,paramterFile)

elif args.mergeLibraries:
    cmd = "python {} {}".format(libmerge, paramterFile)

elif args.preprocess:
    cmd = "python {} '{}' '{}'".format(preprocess, paramterFile, inputFile)

# elif args.masscal:
#     cmd = "python {} {}".format(mass_cal, paramterFile)

elif args.librarySearch:
    cmd = "python {} {} --queue={} --mem={}".format(libsearch, paramterFile, queue, mem)

elif args.librarySearch_w_pp:
    cmd = "python {} {} --queue={} --mem={}".format(libsearch_w_preprocess, paramterFile, queue, mem)

elif args.libraryFilter:
    cmd = "python {} {}".format(libfilter, paramterFile)

elif args.libraryQuan:
    cmd = "python {} {}".format(libquan, paramterFile)
else:
    print ("\n\nPlease select a valid mode.")

try:
    os.system(cmd)
except:
    print (msg(name=None))

