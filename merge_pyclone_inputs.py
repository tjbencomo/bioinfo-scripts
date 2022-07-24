"""
Combine primary and metastatic pyclone input files into a single file for 
pyclone analysis

Usage:
    python merge_pyclone_inputs.py [Primary TSV] [Metastatic TSV] [Output TSV]
"""

import sys
import pandas as pd

def main():
    primaryFp = sys.argv[1]
    metFp = sys.argv[2]
    outFp = sys.argv[3]

    pdf = pd.read_csv(primaryFp, sep = '\t')
    mdf = pd.read_csv(metFp, sep = '\t')
    df = pd.concat([pdf, mdf], ignore_index = True)
    df.to_csv(outFp, index = False, sep = '\t')

if __name__ == '__main__':
    main()
