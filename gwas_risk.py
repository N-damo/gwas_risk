#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import argparse
import os
import subprocess
from argparse import ArgumentParser

def get_args():
    parser=ArgumentParser(description='polygenic prediction')
    parser.add_argument('-vcf','--vcf',help='input your vcf file')
    parser.add_argument('-o','--output',help='working directory')
    args=parser.parse_args()
    return args

def pbs(vcffile, summary_data, id, chinese, output):

    if output.endswith('/') == False:
        output = output + '/'
    assert os.path.exists(output),'error, the {} does not exist,please check it'.format(output)
    with open(id, 'w') as f:
        f.write('#$ -N {}\n'.format(id))
        f.write('#$ -pe smp 3\n')
        f.write('#$ -q all.q\n')
        f.write('#$ -cwd\n')
        f.write('cd {}\n'.format(output))
        f.write('source ~/.bash_profile\n')
        f.write('python /share/data3/lianlin/soft/bin/gwas/gwas_risk2.py -vcf {} -d {} -id {} -sd /share/data2/leon/p3/{}_LDpred_p1.0000e-03.adjusted.txt ;'.format(vcffile, chinese,id,summary_data))

def alljob(jobList):
    with open('poly.bat', 'w') as f:
        for job in jobList:
            f.write('qsub {}\n'.format(job))
    cmd = 'chmod +x poly.bat'
    subprocess.call(cmd, shell = True)


if __name__ == '__main__':
    args = get_args()
    vcffile = args.vcf
    output = args.output
    df = pd.read_table('/share/data2/leon/report/info.txt', header = 0, sep = '\t')
    for i in range(len(df)):
        summary_data = df.loc[i,'dataset']
        id = df.loc[i,'id']
        chinese=df.loc[i,'chinese']
        pbs(vcffile, summary_data, id, chinese,output)
    jobList = df['id']
    alljob(jobList)

    
    