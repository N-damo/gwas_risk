#!/usr/bin/env python
# coding: utf-8


import gzip
import math
import pandas as pd
import argparse
import json
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict
from scipy.stats import norm


def get_args():
    parser=ArgumentParser(description='height prediction')
    parser.add_argument('-vcf','--vcf',help='input your vcf file')
    parser.add_argument('-sd','--summary_data',help='input your gwas summary data')
    parser.add_argument('-d','--disease',help='disease name')
    parser.add_argument('-id','--id',help='disease id')
    args=parser.parse_args()
    return args


def gwas_stat(summary_data):
    gwas = {}
    with open(summary_data,'r') as f:
        next(f)
        for i in f:
            line = i.strip().split('\t')
            rsid = line[1]
            riskallele = line[3]
            beta = line[5]
            pvalue = line[7]
            if float(pvalue) > 5e-2:
                continue
            else:
                freq = line[4]
                gwas[rsid] = [riskallele,beta,pvalue,freq]
    gwas=sorted(gwas.items(),key=lambda x:x[1][2])[0:25]
    gwas=dict(gwas)
    return gwas



def parse_vcf(vcffile, summary_data):
    unvariant = '0/0'
    het = '0/1'
    hom = '1/1'
    gwas = gwas_stat(summary_data)
    genotype = defaultdict(int)
    with gzip.open(vcffile,'r') as f:
        for i in f:
            if i.startswith('##'):
                continue
            elif i.startswith('#CHROM'):
                vcf = i.strip().split('\t')                
                n = len(vcf[9:])
            else: 
                line = i.split('\t')
                rsid = line[2]
                if rsid in gwas:
                    ref = line[3]
                    alt = line[4]
                    riskallele = gwas[rsid][0]
                    freq = float(gwas[rsid][-1])
                    beta = float(gwas[rsid][1])
                    for j in range(n):    
                        risk = 0
                        alleles = line[9 + j].split(':')[0]
                        if alleles == unvariant and ref == riskallele:
                            x = 2
                        elif alleles == het:
                            x = 1
                        elif alleles == hom and alt == riskallele:
                            x = 2
                        else:
                            x = 0
                        #beta*(x-2f)/sqrt(2f(1-f))
                        risk = beta*(x-(2*freq))/math.sqrt(2*freq*(1-freq))
                        genotype[vcf[9+j]] += risk
                      
    return genotype

     
def union_dict(vcffile, summary_data, disease):
    genotype = parse_vcf(vcffile, summary_data)
    df = pd.DataFrame.from_dict(genotype,orient='index')
    df.reset_index(inplace = True)
    df.columns = ['sample',disease]
    df.to_csv(args.id + '.csv',header = True, index = False, encoding='utf_8')
    

if __name__ == '__main__':
    args = get_args()
    vcffile = args.vcf
    summary_data = args.summary_data
    disease = args.disease
    union_dict(vcffile,summary_data, disease)
    


   

   
    
   