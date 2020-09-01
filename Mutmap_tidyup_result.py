#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--idx", help="snp_index result: dir/file")
    parser.add_argument("--eff", help="snpEff result: dir/file")
    parser.add_argument("--out", help="output: dir/file")
    
    args = parser.parse_args()
    
    input1 = str(args.idx)
    input2 = str(args.eff)
    output = str(args.out)
    
    df1 = pd.read_table(f'{input1}', names = ['#CHROM', 'POS', 'VARIANT', 'DEPTH', 'p99','p95','SNP-index'])   
    df2 = pd.read_table(f'{input2}', skiprows=48)
    df3 = df2[['#CHROM','REF','POS','INFO']]
    df4 = df1.merge(df3, how='left', on=['#CHROM','POS'])
    ANN_col1 = df4['INFO'].str.split('ANN=', n = 0, expand = True)
    ANN_col2 = ANN_col1[1].str.split('|', n = 0, expand = True)
    df5 = pd.concat([df4.iloc[:, 0:8],ANN_col2], axis = 1).rename(columns={0: "ALT",
                                                                           1 : "Annotation",
                                                                           2 : "Putative_impact",
                                                                           3 : "Gene Name",
                                                                           4 : "Gene ID",
                                                                           5 : "Feature type",
                                                                           6 : "Feature ID",
                                                                           7 : "Transcript biotype",
                                                                           8 : "Rank/total",
                                                                           9 : "HGVS.c",
                                                                           10 : "HGVS.p",
                                                                           11 : "cDNA_position/cDNA_len",
                                                                           12 : "CDS_position/CDS_len",
                                                                           13 : "Protein_position/Protein_len",
                                                                           14 : "Distance to feature"})
    df5.to_csv(f'{output}',index=False)
