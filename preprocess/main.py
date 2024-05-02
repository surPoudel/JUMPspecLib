import os, sys
from preprocess_functions import main
#main(in_ms2file, dftxt, tmt, ions_test, neutralIonsTest, jump_mod_dict, sta_AA, resultsDirectory)


in1 = sys.argv[1]
in2 = sys.argv[2]
in3 = sys.argv[3]
in4 = sys.argv[4]
in5 = sys.argv[5]
in6 = sys.argv[6]
in7 = sys.argv[7]


main(in1, in2, in3, in4, in5, in6, in7)
