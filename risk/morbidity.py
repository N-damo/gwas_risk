#!/usr/bin/env python3
# coding:utf-8

import pysam
import sys
import os
import json
import math
from collections import defaultdict
import subprocess
import pandas as pd
import re
# import psutil
# import os


# info = psutil.virtual_memory()
# print (u'内存使用：',psutil.Process(os.getpid()).memory_info().rss/(1024**2))
# print (u'总内存：',info.total/(1024**2))
# print (u'内存占比：',info.percent)
# print (u'cpu个数：',psutil.cpu_count())

class Beta(object):
    def __init__(self, vcf, all_summary_data):
        self.vcf = vcf
        self.all_summary_data = all_summary_data
        self.all_dataset = self.gwas_stat()

    def gwas_stat(self):
        #以使用的data作为key，例如
        with open(self.all_summary_data, 'r') as f:
            for i in f:
                line = i.strip()
                all_trait = json.loads(line)
        return all_trait

    def source_dict(self):
        sample_beta = defaultdict(dict)
        vcf = pysam.VariantFile(self.vcf)
        sample_list = list(vcf.header.samples)
        trait_keys = list(self.all_dataset.keys())
        for trait in trait_keys:
            for sample in sample_list:
                sample_beta[trait][sample] = 0
        vcf.close()
        return sample_beta

    def pick_locus(self):
        collect=[]
        for trait in self.all_dataset:
            rsid=self.all_dataset[trait]
            rsid_keys=rsid.keys()
            for id in rsid_keys:
                if id not in collect:
                    collect.append(id)
        with open('gwas-rsid.txt','w') as f:
            f.write('\n'.join(collect))


    def vcf_pick_site(self):
        self.pick_locus()
        cmd='vcftools --gzvcf {vcf} --recode --snps gwas-rsid.txt -c |bgzip > report.vcf.gz ; tabix report.vcf.gz'.format(vcf=self.vcf)
        subprocess.call(cmd,shell=True)


    def parse_variant(self):
        self.vcf_pick_site()
        vcf = pysam.VariantFile('report.vcf.gz')
        sample_beta = self.source_dict()
        # sample_list=list(vcf.header.samples)
        for rec in vcf.fetch():
            rsid = rec.id
            sample_alleles = {
                s.name: list(set(s.alleles)) for s in rec.samples.values()}
            for trait in self.all_dataset:
                summary_data = self.all_dataset[trait]
                if not rsid in summary_data:
                    continue
                else:
                    # print(rsid)
                    #gwas[rsid] = [riskallele, beta, pvalue, freq]
                    riskallele = summary_data[rsid][0]
                    beta = float(summary_data[rsid][1])
                    freq = float(summary_data[rsid][3])
                    for sample in sample_alleles:
                        allele_dup = sample_alleles[sample]
                        if len(allele_dup) == 2:  # double == het genotype
                            x = 1
                        elif allele_dup[0] == riskallele:  # single == ref or alt,if equal to riskallele,x=2
                            x = 2
                        else:
                            x = 0
                        # beta*(x-2f)/sqrt(2f(1-f))
                        glm_beta = beta*(x-(2*freq)) / \
                            math.sqrt(2*freq*(1.0-freq))
                        sample_beta[trait][sample] += glm_beta
        vcf.close()
        return sample_beta

    def infer_gender(self):
        sample_x_hom = defaultdict(int)
        sample_gender = {}
        vcf = pysam.VariantFile(self.vcf)
        total = 0
        for rec in vcf.fetch('X', 2781478, 155701383):
            total += 1
            sample_alleles = {s.name: len(set(s.alleles))
                              for s in rec.samples.values()}
            for sample in sample_alleles:
                #print(sample_alleles[sample])
                if sample_alleles[sample] == 1:
                    sample_x_hom[sample] += 1
        for sample in vcf.header.samples:
            x_hom = sample_x_hom[sample]/total
            if x_hom > 0.98:
                sample_gender[sample] = 0  # male
            else:
                sample_gender[sample] = 1  # female
        vcf.close()
        return sample_gender



    