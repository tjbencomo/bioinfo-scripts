"""
Create input file for pyclone-vi
This script uses a somatic mutations file (MAF file) and copy number alterations
file (info derived from FACETS) and combines the files for use with pyclone-vi.
The end result is a TSV file with the following columns:
    mutation_id
    sample_id
    ref_counts
    alt_counts
    major_cn
    minor_cn
    normal_cn
    tumor_content

This script replaces any LCN_EM values of . with 0, following the
Van Allen Lab's same implementation of FACETS for Terra.
Major copy number is calculated as total copy number (TCN_EM) - lesser copy number (LCN_EM)

The MAF should be from a single sample. If the patient has multiple tumor samples,
run this script on each sample individually and then concatenate. Use the Sample ID flag
to denote whether the sample is primary/metastatic etc

Usage:
    python make_pyclone_input [Patient ID] [Sample ID] [FACETS VCF] [Mutect2 MAF] [FACETS CSV] [Output TSV]
"""
import sys
import os
import gzip
import pandas as pd

def load_maf(fp):
    with open(fp, 'r') as f:
        line = f.readline()
    if line.startswith('#'):
        return pd.read_csv(fp, sep = '\t', skiprows = 1)
    else:
        return pd.read_csv(fp, sep = '\t')


"""
Extract purity and ploidy info from FACETS generated VCF file
"""
def get_purity_ploidy(fp):
    with gzip.open(fp, 'rt') as f:
        line = f.readline()
        purity = None
        ploidy = None
        print("Start searching...")
        while purity is None or ploidy is None:
            if '##purity' in line:
                purity = line.rstrip().replace('##purity=', '')
                print("Found purity!")
            elif '##ploidy' in line:
                ploidy = line.rstrip().replace('##ploidy=', '')
                print("Found ploidy!")
            line = f.readline()
    return purity, ploidy

"""
Given a dataframe with CN segment information, create
a dict mapping segment info to genotype info

Genotype info includes:
    SV type (DUP, DEL, LOH etc)
    Major + minor + normal copy number
"""
def build_cn_lookup(df):
    out = {}
    for index, row in df.iterrows():
        chrom, start, end = row['chrom'], row['start'], row['end']
        start, end = int(start), int(end)
        svtype = row['svtype']
        majorCN, minorCN = row['major_cn'], row['lcn_em']
        if chrom == 'chrY':
            normalCN = 1
        else:
            normalCN = 2
        segmentInfo = (start, end, svtype, majorCN, minorCN, normalCN)
        if chrom in out:
            out[chrom].append(segmentInfo)
        else:
            out[chrom] = [segmentInfo]
    return out

def isInSegment(start, end, segmentInfo):
    # if start == 22766693 and end == 22766693:
    #     print(segmentInfo)
    #     print(f"{start >= segmentInfo[0]}")
    #     print(f"{end <= segmentInfo[1]}")
    return start >= segmentInfo[0] and end <= segmentInfo[1]


"""
Lookup CN segment genotype info for a given mutation that falls
within that segment

position should be formatted as 'chrom:start:end'
"""
def query_cn_segment(position, lookup):
    chrom, start, end = position.split(':')
    start, end = int(start), int(end)
    if chrom not in lookup:
        print(f"Missing chrom {chrom} in lookup table!")
        return None
    # if position == 'chr14:22766693:22766693':
    #         print(lookup[chrom])
    for segment in lookup[chrom]:
        if isInSegment(start, end, segment):
            return segment
    print(f'Could not find a CN segment {position}')
    return None 

"""
Combine mutations from MAF file with CN segments from FACETS
with purity/ploidy information
"""
def merge_alterations(patientId, sampleId, purity, ploidy, mafFp, cnaFp, outFp):
    muts = load_maf(mafFp)
    # muts = muts[muts['Start_Position'] == muts['End_Position']]
    # muts = pd.read_csv(mafFp, sep = '\t', skiprows = 1)
    cnas = pd.read_csv(cnaFp)
    print(muts.head())
    print(cnas.head())

    cnas['lcn_em'] = cnas['lcn_em'].replace(['.'], '0')
    cnas['lcn_em'] = cnas['lcn_em'].astype(int)
    cnas['major_cn'] = cnas['tcn_em'] - cnas['lcn_em']
    cnInfo = build_cn_lookup(cnas)

    muts['mid'] = patientId  + ':' + muts['Chromosome'] + ':' + muts['Start_Position'].astype(str) + ':' + muts['End_Position'].astype(str) + muts['Reference_Allele'] + '>' + muts['Tumor_Seq_Allele2']
    numExcludedMuts = 0
    output = []
    for index, row in muts.iterrows():
        pos = f"{row['Chromosome']}:{row['Start_Position']}:{row['End_Position']}"
        segment = query_cn_segment(pos, cnInfo)
        # if pos == 'chr14:22766693:22766693':
        #     print(segment)
        if segment is None:
            numExcludedMuts += 1
            continue
        mutInfo = (row['mid'], sampleId,  row['t_ref_count'], row['t_alt_count'], segment[3], segment[4], segment[5], purity)
        output.append(mutInfo)
    df = pd.DataFrame(output, columns = ['mutation_id', 'sample_id', 'ref_counts', 
        'alt_counts', 'major_cn', 'minor_cn', 'normal_cn', 'tumour_content'])
    print(f"Not able to find a CN segment for {numExcludedMuts} mutations")
    return df


def main():
    patientId = sys.argv[1]
    sampleId = sys.argv[2]
    vcfFp = sys.argv[3]
    mafFp = sys.argv[4]
    cnaFp = sys.argv[5]
    outFp = sys.argv[6]
    purity, ploidy = get_purity_ploidy(vcfFp)
    print(f"Purity: {purity} Ploidy: {ploidy}")
    df = merge_alterations(patientId, sampleId, purity, ploidy, mafFp, cnaFp, outFp)
    df.to_csv(outFp, sep = '\t', index = False)



if __name__ == "__main__":
    main()
