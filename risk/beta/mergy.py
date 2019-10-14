#!/usr/bin/env python3 
#coding:utf-8

import json
import os
import sys
import re
import numpy as np
"""
合并/share/data2/leon/polygenic/beta和/share/data2/leon/polygenic/sub目录下的json文件
sub：输出trait_locus.json，包含所有数据的位点信息,例如{"1210.assoc.tsv": {"rs1715033": ["C", "6.9479e-02", 0.00255456, "0.0853"], "rs508577": ["A", "6.8299e-02", 0.00549432, "0.0833"]
dir：输出sample_beta_collect.json，包含所有数据的样本beta值，例如{"1960.assoc.tsv": {"E229884": -1.1919094246892652
dir：输出beta_avg_std.json，包含所有数据平均beta值和方差，例如{"1960.assoc.tsv": [-0.05027780171462909, 2.2388760012008797],
注意是在当前目录运行
"""
class Sub_dir(object):
    def __init__(self):
        self.d=self.collect_json()
        json_out(self.d,'trait_locus.json')

    def collect_json(self):
        files=os.listdir('./')
        json_beta=[s for s in files if s.startswith('cut_') and s.endswith('.json')]
        d={}
        for j in json_beta:
            data={}
            with open(j,'r') as f:
                dataset=self.get_dataset(j)
                for i in f:
                    line=i.strip()
                    data[dataset]=json.loads(line)
            d.update(data)
        return d


    def get_dataset(self,j):
        data=os.path.basename(j)
        dataset=re.findall(r'cut_(.*)_LDpred_p1.0000e-03.adjusted.txt.json',data)[0]
        return dataset



class Beta_dir(object):
    def __init__(self):
        self.d=self.collect_json()
        json_out(self.d,'sample_beta_collect.json')
        self.cal_avg_std=self.avg_std()
        json_out(self.cal_avg_std,'beta_avg_std.json')

    def collect_json(self):
        files=os.listdir('./')
        json_beta=[s for s in files if s.startswith('beta') and s.endswith('.json') and s != 'beta_avg_std.json']
        d={}
        for j in json_beta:
            with open(j,'r') as f:
                for i in f:
                    line=i.strip()
                    jd=json.loads(line)
            d.update(jd)
        return d


    def avg_std(self):
        trait=list(self.d.keys())
        cal_avg_std={}
        for t in trait:
            values=self.d[t]
            beta=list(values.values())
            avg=np.mean(beta)
            std=np.std(beta,ddof=1)
            cal_avg_std[t]=[avg,std]
        return cal_avg_std



def json_out(d,out_name):
    with open(out_name,'w') as f:
        json.dump(d,f)
        f.write('\n')


if __name__ == '__main__':
    directory=sys.argv[1]
    if directory == 'sub':
        if os.path.basename(directory) != 'sub':
            print('are you sure in sub?')
            sys.exit()
        Sub_dir()
    elif directory == 'beta':
        if os.path.basename(directory) != 'beta':
            print('are you sure in beta?')
            sys.exit()
        Beta_dir()
    else:
        print('sub_dir or beta_dir,not else')
        sys.exit()
