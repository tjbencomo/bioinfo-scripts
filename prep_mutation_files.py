"""
Given a set of primary and metastatic tumor MAF files, figure out which mutations
are exclusive to each sample. These exclusive sites need to be present in both samples,
but with a alt count of 0 in the sample without the SNV but still requires the reference 
read count in the non-mutant sample. This script gets the number of reference reads in the non-mutant
sample, adds the missing SNV to the non-mutant MAF file with an alt_count of 0, and saves the MAF
file for use with make_pyclone_input. 

This process is repeated for both MAF files

Usage:
    python prep_mutation_files.py [Primary MAF] [Primary Tumor Bam] \
            [Metastatic MAF] [Metastatic Tumor Bam] \
            [Output directory]
"""

import os
import sys
import pandas as pd
import pysam as ps
import pysamstats as pst
import subprocess as sb

"""
pdf - primary sample MAF dataframe
mdf - metastatic sample MAF dataframe
"""
def find_exclusive_mutations(pdf, mdf):
    pdf['mid'] = pdf['Chromosome'] + ':' + pdf['Start_Position'].astype(str) + ':' + pdf['End_Position'].astype(str) + ':' + pdf['Reference_Allele'] + '>' + pdf['Tumor_Seq_Allele2']
    
    mdf['mid'] = mdf['Chromosome'] + ':' + mdf['Start_Position'].astype(str) + ':' + mdf['End_Position'].astype(str) + ':' +  mdf['Reference_Allele'] + '>' + mdf['Tumor_Seq_Allele2']

    pMutIds =  set(pdf['mid'].tolist())
    mMutIds = set(mdf['mid'].tolist())
    pSpecificIds = pMutIds - mMutIds
    mSpecificIds = mMutIds - pMutIds

    return pSpecificIds, mSpecificIds

"""
mutdf - mutation MAF dataframe for a sample
missingdf - list of mutation IDs not present in the mutdf sample

This script uses pysamstats to get the number of reference reads for each
mutation not present in the sample
"""
def add_absent_mutations(mutdf, absentMuts, bamfp):
    bam = ps.AlignmentFile(bamfp, 'rb')
    records = []
    N = len(absentMuts)
    skipCounter = 0
    for i, mut in enumerate(absentMuts):
        if i % 1000 == 0:
            print(f"Finished querying {i} mutations")
        # print(f"Running pysamstats for {mut}")
        chrom, start, end, alleles  = mut.split(':')
        ref, alt = alleles.split('>')
        info = pst.stat_coverage(bam, chrom = chrom, start = int(start), end = int(end) + 1, 
                truncate = True, one_based = True)
        info = list(info)
        if len(info) == 0:
            # Very rarely no results - skip mutation
            skipCounter += 1
            continue
        info = info[0]
        d = {
            'Chromosome' : chrom,
            'Start_Position' : start,
            'End_Position' : end,
            'Reference_Allele' : ref,
            'Tumor_Seq_Allele2' : alt,
            't_ref_count' : info['reads_all'],
            't_alt_count' : 0
        }
        records.append(d)
    print(f"Skipped {skipCounter} mutations while querying from BAM")
    newdf = pd.DataFrame(records)
    return pd.concat([mutdf, newdf], axis = 0, ignore_index = True)

def load_maf(fp):
    with open(fp, 'r') as f:
        line = f.readline()
    if line.startswith('#'):
        return pd.read_csv(fp, sep = '\t', skiprows = 1)
    else:
        return pd.read_csv(fp, sep = '\t')

def main():
    pMafFp = sys.argv[1]
    pBamFp = sys.argv[2]
    mMafFp = sys.argv[3]
    mBamFp = sys.argv[4]
    outDir = sys.argv[5]
   
    print("Loading MAF files")
    primaryMaf = load_maf(pMafFp)
    metMaf = load_maf(mMafFp)
    # primaryMaf = pd.read_csv(pMafFp, sep = '\t')
    # metMaf = pd.read_csv(mMafFp, sep = '\t')
    print("Finding exclusive mutations")
    pSpecific, mSpecific = find_exclusive_mutations(primaryMaf, metMaf)
    print(f"There are {len(pSpecific)} primary-specific mutations and {len(mSpecific)} metastatsis-specific mutations")

    print("Adding absent mutations to primary MAF")
    primaryMaf = add_absent_mutations(primaryMaf, mSpecific, pBamFp)
    print("Adding absent mutations to metastatic MAF")
    metMaf = add_absent_mutations(metMaf, pSpecific, mBamFp)
    
    pOutFn = pMafFp.replace('.maf', '')
    pOutFn = f"{pOutFn}.primary.with_absent_mutations.maf"
    pOutFp = os.path.join(outDir, os.path.basename(pOutFn))
    mOutFn = mMafFp.replace('.maf', '')
    mOutFn = f"{mOutFn}.metastatis.with_absent_mutations.maf"
    mOutFp = os.path.join(outDir, os.path.basename(mOutFn))
    print("Saving updated MAF files")
    primaryMaf.to_csv(pOutFp, index = False, sep = '\t')
    metMaf.to_csv(mOutFp, index = False, sep = '\t')
    

if __name__ == '__main__':
    main()
