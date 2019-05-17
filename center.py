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


def ukbiobank(sample,hg38):
    'variant rsid    nCompleteSamples        AC      ytx     beta    se      tstat   pval'
    
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    id = {}
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[1]           
            if rsid not in id:
                id[rsid] = 0
                if rsid in hg38:
                    chr = hg38[rsid].split(':')[0]
                    position = hg38[rsid].split(':')[1]
                    ref_a2 = i[0].split(':')[-2]
                    alt_a1 = i[0].split(':')[-1]
                    beta = i[-4]
                    if is_number(beta):                   
                        odd_ratio = str(math.e**float(beta))
                        pvalue = i[-1]
                        if 0<float(pvalue)<1:
                            f.write(chr + '\t')
                            f.write(rsid + '\t')
                            f.write(alt_a1 + '\t')
                            f.write(ref_a2 + '\t')
                            f.write(position + '\t')
                            f.write(odd_ratio + '\t')
                            f.write(pvalue + '\n')
                    else:
                        print 'dirty data'
                        cmd1 = 'echo {} >>dirty_data'.format(output)
                        
                        subprocess.call(cmd1, shell=True)
                    
                   
               
def magenetic(sample,hg38):
    'chromosome position ID EA NEA eaf beta se p-value n_studies n_samples'
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            gwas.append(i.strip().split(' '))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[2]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4]
                alt_a1 = i[3]
                beta = i[6]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[8]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
                    
                    subprocess.call(cmd1, shell=True)
                
                    
def dbgap(sample,hg38):
    """SNP ID  P-value Chr ID  Chr Position    Submitted SNP ID        ss2rs   rs2genome       Allele1 Allele2 pHWE (case)     pHWE (control)  Ca
ll rate (case)        Call rate (control)     Odds ratio      CI low  CI high"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[7]
                alt_a1 = i[8]
                odd_ratio = i[-3]
                if is_number(odd_ratio):                   
                    pvalue = i[1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
                  
                    subprocess.call(cmd1, shell=True)
                    
                   
                    
def eagle(sample,hg38):
    """rsID    chromosome      position        reference_allele        other_allele    eaf     European_N      beta    se      p.value i2      q_
p.value       AllEthnicities_N        log10BF"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4]
                alt_a1 = i[3]
                beta = i[7]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[9]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
                    
                    subprocess.call(cmd1, shell=True)
                   
                    


def egg(sample,hg38):
    """CHR     BP      RSID    EA      NEA     BETA    SE      P       N"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[2]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[5]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[7]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)
                   
def gefos(sample,hg38):
    """chromosome      position        rs_number       reference_allele        other_allele    eaf     beta    se      beta_95L        beta_95U        
    z       p-value _-log10_p-value q_statistic     q_p-value       i2      n_studies       n_samples       effects"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[2]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4]
                alt_a1 = i[3]
                beta = i[6]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[11]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)



def giant(sample,hg38):
    """MarkerName Allele1 Allele2 Allele1_Freq_HapMapCEU beta se P N_cases N_controls"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split(' '))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[2].upper()
                alt_a1 = i[1].upper()
                beta = i[4]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[6]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def gis(sample,hg38):
    """SNP CHR BP A1 A2 FREQ_A1 EFFECT_A1 SE P"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split(' '))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4]
                alt_a1 = i[3]
                beta = i[-3]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:

                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def glgc(sample,hg38):
    """SNP_hg18        SNP_hg19        rsid    A1      A2      beta    se      N       P-value Freq.A1.1000G.EUR"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[2]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[5]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-2]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def gugc(sample,hg38):
    """MarkerName,n_total,A1,A2,beta,se,p_gc"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split(','))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[3].upper()
                alt_a1 = i[2].upper()
                beta = i[4]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def ibdgc(sample,hg38):
    """SNP Chr Pos A1_effect A2_other Mantra_n_studies Mantra_log10BF Mantra_n_samples Mantra_dir 
    beta_EUR se_EUR P_EUR beta_EAS se_EAS pval_EAS beta_IND se_IND pval_IND beta_IRA se_IRA pval_IRA EAF_EUR EAF_EAS EAF_IND EAF_IRA"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split(' '))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[12]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[14]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[14]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)


def kilpelainen(sample,hg38):
    """SNP     Effect  StdErr  P       N       Effect_all      Other_all"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[-1].upper()
                alt_a1 = i[-2].upper()
                beta = i[1]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[3]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[3]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)


def lanegroup(sample,hg38):
    """SNP CHR BP A1 A2 MAF BETA SE P N"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split(' '))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[6]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-2]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-2]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
               

def magic(sample,hg38):
    """snp     effect_allele   other_allele    maf     effect  stderr  pvalue"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[2].upper()
                alt_a1 = i[1].upper()
                beta = i[-3]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def na(sample,hg38):
    """CHR MARKERNAME POS EFFECT_ALLELE OTHER_ALLELE BETA STDERR P DIRECTION HET_P"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split(' '))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[1]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[5]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-3]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-3]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print beta
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def ssgac(sample,hg38):
    """MarkerName      CHR     POS     A1      A2      EAF     Beta    SE      Pval"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split(' '))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[6]
                if is_number(beta):                   
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print beta
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)


def pgc(sample,hg38):
    """CHR     SNP     BP      A1      A2      INFO    OR      SE      P"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[1]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[6]
                if is_number(beta):                   
                    odd_ratio = str(beta)
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-1]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print beta
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def pgc_large(sample,hg38):
    """CHR     SNP     BP      A1      A2      FRQ_A_5305      FRQ_U_5305      INFO    OR      SE      P       ngt"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[1]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[8]
                if is_number(beta):                   
                    odd_ratio = str(beta)
                    pvalue = i[-2]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                elif 'e' in beta:
                    odd_ratio = str(math.e**float(beta))
                    pvalue = i[-2]
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print beta
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)

def pgc_small(sample,hg38):
    """snpid hg18chr bp a1 a2 or se pval info ngt CEUaf"""
    gwas = []
    with gzip.open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/center/basic2/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[4].upper()
                alt_a1 = i[3].upper()
                beta = i[5]
                if is_number(beta) and float(beta) > 0:                   
                    odd_ratio = str(beta)
                    pvalue = i[7]
                    if 0<float(pvalue)<1:
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
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print beta
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)


def cancer1(sample,hg38):
    """id      a1      a2      or      p"""
    gwas = []
    with open(sample, 'r') as f:
        for i in f:
            if i.startswith('#'):
                continue
            else:
                gwas.append(i.strip().split('\t'))
    
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    print output
    with open('/share/data2/leon/gwas/catalog/basic3/' + output,'w') as f:
        for i in gwas:
            rsid = i[0]
            if rsid in hg38:
                chr = hg38[rsid].split(':')[0]
                position = hg38[rsid].split(':')[1]
                ref_a2 = i[2].upper()
                alt_a1 = i[1].upper()
                beta = i[3]
                if is_number(beta) and float(beta) > 0:                   
                    odd_ratio = str(beta)
                    pvalue = i[4]
                    if 0<float(pvalue)<1:
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
                    if 0<float(pvalue)<1:
                        f.write(chr + '\t')
                        f.write(rsid + '\t')
                        f.write(alt_a1 + '\t')
                        f.write(ref_a2 + '\t')
                        f.write(position + '\t')
                        f.write(odd_ratio + '\t')
                        f.write(pvalue + '\n')
                else:
                    print beta
                    print 'dirty data'
                    cmd1 = 'echo {} >>dirty_data'.format(output)
             
                    subprocess.call(cmd1, shell=True)


def run(data_format):
    if data_format == 'ukbiobank':
        ukbiobank(sample,hg38)
    elif data_format == 'magenetic':
        magenetic(sample,hg38)
    elif data_format == 'dbgap':
        dbgap(sample,hg38)
    elif data_format == 'eagle':
        eagle(sample,hg38)
    elif data_format == 'egg':
        egg(sample,hg38)
    elif data_format == 'gefos':
        gefos(sample,hg38)
    elif data_format == 'giant':
        giant(sample,hg38)
    elif data_format == 'gis':
        gis(sample,hg38)
    elif data_format == 'glgc':
        glgc(sample,hg38)
    elif data_format == 'gugc':
        gugc(sample,hg38)
    elif data_format == 'ibdgc':
        ibdgc(sample,hg38)
    elif data_format == 'kilpelainen':
        kilpelainen(sample,hg38)
    elif data_format == 'lanegroup':
        lanegroup(sample,hg38)
    elif data_format == 'magic':
        magic(sample,hg38)
    elif data_format == 'na':
        na(sample,hg38)
    elif data_format == 'ssgac':
        ssgac(sample,hg38)
    elif data_format == 'pgc':
        pgc(sample,hg38)
    elif data_format == 'pgc_large':
        pgc_large(sample,hg38)
    elif data_format == 'pgc_small':
        pgc_small(sample,hg38)
    elif data_format == 'cancer1':
        cancer1(sample,hg38)
    else:
        print 'select right format'
        sys.exit()

if __name__ == '__main__':
    args = get_args()
    sample = args.input
    hg38 = hg38dict()
    data_format = args.format
    run(data_format)