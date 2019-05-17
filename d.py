#!/share/data3/lianlin/soft/python/bin/python2.7
import gzip
import sys
import pandas as pd
import argparse
import subprocess
import os
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description='gwas summary data adjusted by ldpred')
    parser.add_argument('-i', '--input', help='input your file')
    args = parser.parse_args()
    return args


def orig_summary(sample):
    gwas = {}
    with open(sample, 'r') as f:
        for i in f:
            h = i.strip().split('\t')
            rsid = h[1]
            pvalue = h[-1]
            gwas[rsid] = pvalue
    return gwas


def freq():
    chinese = {}
    with open('/share/data2/leon/chinese', 'r') as f:
        for i in f:
            h = i.strip().split('\t')
            rsid = h[2]
            frequency = h[-1]
            ref = h[3]
            alt = h[4]
            chinese[rsid] = ':'.join([ref, alt, frequency])
    return chinese


def ldpred(output):
    with open('/share/data2/leon/gwas/ldpred/ldpred/' + output + '_LDpred_p1.0000e-03.txt', 'r') as f:
        inf = f.readlines()
    return inf


def final_result(inf, gwas, chinese, output):
    with open('/share/data2/leon/p3all/'+ output +'_LDpred_p1.0000e-03.adjusted.txt', 'w') as f:
        f.write(
            'Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbJ\tbJ_se\tpJ\tLD_r\n')
        for i in inf:
            line = i.split()
            rsid = line[2]
            beta = line[-1]
            a1 = line[3]
            if rsid in gwas:
                if a1 == chinese[rsid].split(':')[0]:
                    f.write('-' + '\t')
                    f.write(rsid + '\t')
                    f.write('-' + '\t')
                    f.write(a1 + '\t')
                    f.write(
                        str(1-float(chinese[rsid].split(':')[-1])) + '\t')
                    f.write(beta + '\t')
                    f.write('-' + '\t')
                    f.write(gwas[rsid] + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\n')
                elif a1 == chinese[rsid].split(':')[1]:
                    f.write('-' + '\t')
                    f.write(rsid + '\t')
                    f.write('-' + '\t')
                    f.write(a1 + '\t')
                    f.write(chinese[rsid].split(':')[-1] + '\t')
                    f.write(beta + '\t')
                    f.write('-' + '\t')
                    f.write(gwas[rsid] + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\t')
                    f.write('-' + '\n')




if __name__ == '__main__':
    args = get_args()
    sample = args.input
    output = sample.split('/')[-1]
    gwas = orig_summary(sample)
    chinese = freq()
    inf = ldpred(output)
    final_result(inf, gwas, chinese, output)


