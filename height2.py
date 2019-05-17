#!/usr/bin/env python
# coding: utf-8


import gzip
import math
import pandas as pd
import argparse
import json
from argparse import ArgumentParser
from collections import defaultdict
import chardet


def get_args():
    parser = ArgumentParser(description='height prediction')
    parser.add_argument('-vcf', '--vcf', help='input your vcf file')
    parser.add_argument('-sd', '--summary_data',
                        help='input your gwas summary data')
    parser.add_argument('-height', '--height', help='input your height info')
    parser.add_argument(
        '-p', '--pvalue', help='select pvalue threshold', default='1e-7')
    parser.add_argument('-result', '--result', help='result prefix')
    args = parser.parse_args()
    return args


class beta(object):

    def __init__(self, vcf, height, summary_data):

        self.vcf = vcf
        self.height = self.sample_info(height)
        self.gwas = self.gwas_stat(summary_data)
        self.parse_vcf()

    def gwas_stat(self, summary_data):
        gwas = {}
        with open(summary_data, 'r') as f:
            next(f)
            for i in f:
                line = i.strip().split('\t')
                rsid = line[1]
                riskallele = line[3]
                beta = line[5]
                pvalue = float(line[7])
                freq = line[4]
                gwas[rsid] = [riskallele, beta, pvalue, freq]
        gwas = sorted(gwas.items(), key=lambda x: abs(
            float(x[1][1])), reverse=True)[0:16]
        gwas = dict(gwas)
        for id in gwas:
            if gwas[id][-2] > 0.05:
                continue
            else:
                gwas.pop('id')
        with open('height.json', 'w') as f:
            json.dump(gwas, f)
        print 'ok'
        return gwas

    def sample_info(self, height):
        _ = pd.read_csv(height, header=0)
        return _

    def parse_vcf(self):
        unvariant = '0/0'
        het = '0/1'
        hom = '1/1'
        genotype = defaultdict(int)
        with gzip.open(self.vcf, 'r') as f:
            for i in f:
                if i.startswith('##'):
                    continue
                elif i.startswith('#CHROM'):
                    vcf = i.strip().split('\t')
                    n = len(vcf[9:])
                else:
                    line = i.split('\t')
                    rsid = line[2]
                    if rsid in self.gwas:
                        ref = line[3]
                        alt = line[4]
                        riskallele = self.gwas[rsid][0]
                        freq = float(self.gwas[rsid][-1])
                        pvalue = float(self.gwas[rsid][-2])
                        beta = float(self.gwas[rsid][1])
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
                                math.sqrt(2*freq*(1-freq))
                            genotype[vcf[9+j]] += risk

        df = pd.DataFrame.from_dict(genotype, orient='index')
        df.columns = ['risk']
        df.reset_index(inplace=True)
        info = pd.merge(self.height, df, on='index')
        info['rank'] = ''
        info['predict'] = ''
        info['y'] = ''
        for i in range(len(info)):
            if info.loc[i, '籍贯'].startswith(('福建', '广东', '江西', '浙江', '广西', '湖南', '湖北', '贵州', '台湾', '云南', '海南')):
                if info.loc[i, '性别'] == '男':
                    info.loc[i, 'rank'] = '南方男'
                else:
                    info.loc[i, 'rank'] = '南方女'
            else:
                if info.loc[i, '性别'] == '男':
                    info.loc[i, 'rank'] = '北方男'
                else:
                    info.loc[i, 'rank'] = '北方女'
            h = get_value()
            h.height = info.loc[i, '身高']
            h.rank = info.loc[i, 'rank']
            h.risk = info.loc[i, 'risk']
            info.loc[i, 'y'] = h.y()
            info.loc[i, 'predict'] = h.predict()
        info.to_csv(args.result + '.csv', header=True,
                    index=False, encoding='utf_8_sig')


class get_value(object):
    avg1_map = {"北方男": 171.476667, "北方女": 163.113333,
                "南方男": 168.674167, "南方女": 160.228333}
    sd1_map = {"北方男": 5.7, "北方女": 5.2, "南方男": 5.7, "南方女": 5.2}

    def y(self):
        return (self.height-self.avg1_map[self.rank])/self.sd1_map[self.rank]

    def predict(self):
        return self.risk*self.sd1_map[self.rank] + self.avg1_map[self.rank]


if __name__ == '__main__':
    args = get_args()
    beta(args.vcf, args.height, args.summary_data)
