#!/usr/bin/python2.7
import gzip
import sys
import math
import argparse
import subprocess
import re
from argparse import ArgumentParser
def get_args():
    parser=ArgumentParser(description='gwas summary data format')
    parser.add_argument('-i','--input',help='input your file')
    parser.add_argument('-format','--format',help='gwas summary data format')
    args=parser.parse_args()
    return args

def is_number(num):
    pattern = re.compile(r'^[-+]?[-0-9]\d*\.\d*|[-+]?\.?[0-9]\d*$')
    result = pattern.match(num)
    if result:
        return True
    else:
        return False

def hg38dict():
    l = []
    with open('/share/data2/leon/hg38','r') as f:
        for i in f:
            l.append(i.strip().split('\t'))

    hg38 = {}
    for i in range(len(l)):
        hg38[l[i][-1]]=':'.join(l[i][:2])
    return hg38

def catalog1(sample,hg38):
    """CHROM  POS     ID      REF     ALT     ALT_FREQ        TEST    OBS_CT  OR      SE      T_STAT  P"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[2]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[3].upper()
                alt_a1 = i[4].upper()
                odd_ratio = i[8]
                if is_number(odd_ratio):                   
                    odd_ratio = str(odd_ratio)
                    pvalue = i[11]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')
                

def catalog2(sample,hg38):
    """regional.analysis       CHR     BP      rsid    A1      A2      A1.AF   P       BETA    SE"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[3]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[5].upper()
                alt_a1 = i[4].upper()
                try:
                    if len(i) == 10:
                        beta = i[8]
                        if is_number(beta):                   
                            odd_ratio = str(math.e**float(beta))
                            pvalue = i[7]
                            if float(pvalue) == 0:
                                pvalue = str(1e-300)
                            f.write(chr + '\t')
                            f.write(rsid + '\t')
                            f.write(alt_a1 + '\t')
                            f.write(ref_a2 + '\t')
                            f.write(position + '\t')
                            f.write(odd_ratio + '\t')
                            f.write(pvalue + '\n')
                        elif 'e' in beta:
                            odd_ratio = str(math.e**float(beta))
                            pvalue = i[7]
                            if float(pvalue) == 0:
                                pvalue = str(1e-300)
                            f.write(chr + '\t')
                            f.write(rsid + '\t')
                            f.write(alt_a1 + '\t')
                            f.write(ref_a2 + '\t')
                            f.write(position + '\t')
                            f.write(odd_ratio + '\t')
                            f.write(pvalue + '\n')
                except:
                    print len(i),i

def catalog3(sample,hg38):
    """Chr Pos Allele1 Allele2 rsID Effect P Direction HetISq N"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[4]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[3].upper()
                alt_a1 = i[2].upper()
                beta = i[5]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[6]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[6]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')


def catalog4(sample,hg38):
    """CHR         SNP         BP   A1      F_A      F_U   A2       
     CHISQ            P           OR           SE          L95          U95"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[1]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[6].upper()
                alt_a1 = i[3].upper()
                odd_ratio = i[9]
                if is_number(odd_ratio) and float(odd_ratio) > 0:                   
                    odd_ratio = str(odd_ratio)
                    pvalue = i[8]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')
                elif 'e' in odd_ratio:
                    odd_ratio = str(odd_ratio)
                    pvalue = i[8]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')

def freq():
    chinese = {}
    with open('/share/data2/leon/chinese', 'r') as f:
        for i in f:
            h = i.strip().split('\t')
            chr = h[0]
            position = h[1]
            rsid = h[2]
            frequency = h[-1]
            ref = h[3]
            alt = h[4]
            chinese[rsid] = ':'.join([chr,position,ref,alt,frequency])
    return chinese

def catalog5(sample,chinese = freq()):
    """28806749-GCST004833-EFO_0001061-Build37.f.tsv"""
    """chromosome      variant_id      base_pair_location      effect_allele   test    nmiss   odds_ratio      standard_error  ci_lower        
    ci_upper        stat    p_value beta    other_allele    effect_allele_frequency"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[1]
            if rsid in chinese:
                chr = chinese[rsid].split(':')[0]
                position = chinese[rsid].split(':')[1]
                alt_a1 = i[3].upper()
                if alt_a1 == chinese[rsid].split(':')[2]:
                    ref_a2 = chinese[rsid].split(':')[3]
                elif alt_a1 == chinese[rsid].split(':')[3]:
                    ref_a2 = chinese[rsid].split(':')[2]
                if len(i) == 15:
                    odd_ratio = i[6]
                    if is_number(odd_ratio):                   
                        odd_ratio = str(odd_ratio)
                        pvalue = i[11]
                        if float(pvalue) == 0:
                            pvalue = str(1e-300)
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                    elif 'e' in odd_ratio:
                        odd_ratio = str(odd_ratio)
                        pvalue = i[11]
                        if float(pvalue) == 0:
                            pvalue = str(1e-300)
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')


def catalog6(sample,hg38):
    """variant_id      chromosome      base_pair_location      effect_allele   other_allele    effect_allele_frequency R2_oncoarray    
    odds_ratio      standard_error  chi2    p_value beta    ci_lower        ci_upper"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[7]
                if is_number(beta) and beta != 'NA':                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[10]
                    if 0 < float(pvalue) <1:
                        #pvalue = str(1e-300)
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[10]
                    if 0 < float(pvalue) <1:
                        #pvalue = str(1e-300)
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')


def catalog7(sample,hg38):
    """chromosome      base_pair_location      variant_id      other_allele    effect_allele   p_value beta    standard_error  
    odds_ratio      ci_lower        ci_upper        effect_allele_frequency"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[2]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[3].upper()
                alt_a1 = i[4].upper()
                beta = i[6]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[5]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[5]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')


def catalog8(sample,hg38):
    """chr     snpid   pos     a1      a2      pvalue  or      se"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split())
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[2]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[3].upper()
                alt_a1 = i[4].upper()
                beta = i[6]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[5]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[5]
                    if float(pvalue) == 0:
                        pvalue = str(1e-300)
                    f.write(chr + '\t')
                    f.write(rsid + '\t')
                    f.write(alt_a1 + '\t')
                    f.write(ref_a2 + '\t')
                    f.write(position + '\t')
                    f.write(odd_ratio + '\t')
                    f.write(pvalue + '\n')


def run(data_format):
    if data_format == 'catalog1':
        catalog1(sample,hg38)
    elif data_format == 'catalog2':
        catalog2(sample,hg38)
    elif data_format == 'catalog3':
        catalog3(sample,hg38)
    elif data_format == 'catalog4':
        catalog4(sample,hg38)
    elif data_format == 'catalog5':
        catalog5(sample,chinese = freq())
    elif data_format == 'catalog6':
        catalog6(sample,hg38)
    elif data_format == 'catalog7':
        catalog7(sample,hg38)



    else:
        print 'select right format'
        sys.exit()

if __name__ == '__main__':
    args = get_args()
    sample = args.input
    hg38 = hg38dict()
    data_format = args.format
    run(data_format)