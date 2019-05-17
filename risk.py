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
    parser = ArgumentParser(description='height prediction')
    parser.add_argument('-vcf', '--vcf', help='input your vcf file')
    parser.add_argument('-sd', '--summary_data',
                        help='input your gwas summary data')
    parser.add_argument(
        '-p', '--pvalue', help='select pvalue threshold', default='5e-6')
    parser.add_argument('-d', '--disease', help='disease name')
    args = parser.parse_args()
    return args


def gwas_stat(summary_data):
    gwas = {}
    with open(summary_data, 'r') as f:
        next(f)
        for i in f:
            line = i.strip().split('\t')
            rsid = line[1]
            riskallele = line[3]
            beta = line[5]
            pvalue = float(line[7])
            if pvalue > 5e-2:
                continue
            else:
                freq = line[4]
                gwas[rsid] = [riskallele, beta, pvalue, freq]
    gwas=sorted(gwas.items(),key=lambda x:abs(float(x[1][1])),reverse=True)[0:100]
    gwas=dict(gwas)
    value=max(gwas,key=lambda x:abs(float(gwas[x][1])))
    print value,gwas[value]
    with open('gwas.json','w') as f:
        json.dump(gwas,f)
    print 'ok'
    return gwas


def parse_vcf(vcffile, p, summary_data):
    unvariant = '0/0'
    het = '0/1'
    hom = '1/1'
    gwas = gwas_stat(summary_data)
    genotype = defaultdict(int)
    with gzip.open(vcffile, 'r') as f:
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
                        # beta*(x-2f)/sqrt(2f(1-f))
                        risk = beta*(x-(2*freq)) / \
                            math.sqrt(2*freq*(1.0-freq))
                        genotype[vcf[9+j]] += risk

    return genotype


def risk_area(genotype):
    risk = {}
    beta = genotype
    values = beta.values()
    avg = np.mean(values)
    std = np.std(values, ddof=1)
    t = norm.ppf(1-(1e-5))#人群发病率为1e-5
    print avg, std
    for sample in beta:
        x = (beta[sample] - avg)/std
        _ = 1-norm.cdf(t-x)
        risk[sample] = _
    with open('{disease}.json'.format(disease=args.disease), 'w') as f:
        json.dump(risk, f)


if __name__ == "__main__":
    args = get_args()
    vcffile = args.vcf
    p = args.pvalue
    summary_data = args.summary_data
    genotype = parse_vcf(vcffile, p, summary_data)
    risk_area(genotype)
