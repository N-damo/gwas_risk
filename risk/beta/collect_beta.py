#!/usr/bin/env python3
#coding:utf-8

import pysam
import math
import re
import sys
import os
import json
from collections import defaultdict

"""
这个脚本是收集样本的beta值，汇总后计算avg和std
输入是ld adjusted后的summary_data，并且是summary_data.py提取过位点后的json文件，带有cut_前缀(数据存放在/share/data2/leon/polygenic/sub)
输出文件在/share/data2/leon/polygenic/beta
"""


class Beta(object):
    def __init__(self, vcf, summary_data):
        self.vcf = vcf
        self.summary_data = summary_data
        self.gwas_stat()

    def gwas_stat(self):
        self.trait = ''
        with open(self.summary_data, 'r') as f:
            for i in f:
                line = i.strip()
                self.trait = json.loads(line)
        return self.trait

    def parse_variant(self):
        sample_beta = defaultdict(int)
        vcf = pysam.VariantFile(self.vcf)
        for rec in vcf.fetch():
            rsid = rec.id
            if not rsid in self.trait:
                continue
            else:
                #gwas[rsid] = [riskallele, beta, pvalue, freq]
                riskallele = self.trait[rsid][0]
                beta = float(self.trait[rsid][1])
                freq = float(self.trait[rsid][3])
                sample_alleles = {
                    s.name: s.alleles for s in rec.samples.values()}
                for sample in sample_alleles:
                    glm_beta = 0
                    allele_dup = set(sample_alleles[sample])
                    if len(allele_dup) == 2:  # double == het genotype
                        x = 1
                    elif allele_dup == riskallele:  # single == ref or alt,if equal to riskallele,x=2
                        x = 2
                    else:
                        x = 0
                    # beta*(x-2f)/sqrt(2f(1-f))
                    glm_beta = beta*(x-(2*freq)) / math.sqrt(2*freq*(1.0-freq))
                    sample_beta[sample] += glm_beta
        vcf.close()
        return sample_beta


    def get_dataset_name(self):
        # """Summary_statistics_MAGNETIC_XL.VLDL.P.txt_LDpred_p1.0000e-03.adjusted.txt"""
        trait_name = os.path.basename(self.summary_data)
        dataset = re.findall(r'^cut_(.*)\.json$', trait_name)[0]
        dataset_name = re.sub('_LDpred_p1.0000e-03.adjusted.txt', '', dataset)
        return dataset_name



def json_out(trait_beta_collect,out_name):
    with open("/share/data2/leon/polygenic/beta/"+out_name,'w') as f:
        json.dump(trait_beta_collect,f)
        f.write('\n')



if __name__ == '__main__':
    trait_beta_collect=defaultdict(dict)
    vcf=sys.argv[1]
    summary_data=sys.argv[2]
    result=Beta(vcf,summary_data)
    sample_beta=result.parse_variant()
    dataset_name=result.get_dataset_name()
    out_name='beta_'+dataset_name+'.json'
    if os.path.exists(out_name):
        with open("/share/data2/leon/polygenic/beta/"+out_name,'r') as f:
            for i in f:
                line=i.strip()
                trait_beta_collect=json.loads(line)
        for sample in sample_beta:
            if sample in trait_beta_collect[dataset_name]:
                continue
            else:
                trait_beta_collect[dataset_name][sample]=sample_beta[sample]
    else:
        for sample in sample_beta:
            trait_beta_collect[dataset_name][sample]=sample_beta[sample]
    json_out(trait_beta_collect,out_name)

    
