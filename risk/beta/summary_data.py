#!/usr/bin/env python3
#coding:utf-8

import sys
import os
import json

"""
提取summary_data有效位点信息，summary_data数据在/share/data2/leon/p3all
输出的json文件在/share/data2/leon/polygenic/sub
"""

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
    print(value,gwas[value])
    out_name='cut_'+os.path.basename(summary_data)+'.json'
    with open('/share/data2/leon/polygenic/sub/'+out_name,'w') as f:
        json.dump(gwas,f)

if __name__ == '__main__':
    summary_data=sys.argv[1]
    gwas_stat(summary_data)