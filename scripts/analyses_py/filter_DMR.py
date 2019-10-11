import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
dmr_file=sys.argv[1]+"_rms_results_collapsed.tsv"
output=sys.argv[1]+"_filtered_rms_results_collapsed.tsv"
min_dms=4
min_mC_diff=0.2

#filter DMRs
functions.filter_dmr(dmr_file,output=output,min_dms=min_dms,min_mC_diff=min_mC_diff)

