#!/share/data3/lianlin/soft/python/bin/python2.7
import gzip
import sys
import pandas as pd
import argparse
import subprocess
import os
from argparse import ArgumentParser


def get_args():
    parser=ArgumentParser(description='gwas summary data adjusted by ldpred')
    parser.add_argument('-i','--input',help='input your file')
    parser.add_argument('-s','--source',help='study or share center')
    parser.add_argument('-radius','--radius',help='LD radius',default = '2000')
    args=parser.parse_args()
    return args




# def hsdf5(sample,output,samplesize):
#     if os.path.exists(output):
#         print '{} exists'.format(output)
#         cmd = 'rm {}*'.format(output)
#         print cmd
#         subprocess.call(cmd, shell = True)
        
#     cmd = "coord --gf=/share/data2/leon/gwas/ldpred/LD_reference --ssf={} --N={} --out={} --ssf_format=BASIC ;".format(sample,samplesize,output)
#     print cmd
#     subprocess.call(cmd, shell = True)
    


def adjusted(output,samplesize,radius):
    cmd = "ldpred --coord={} --PS=0.01,0.0001,0.00001 --ld_radius={} --local_ld_file_prefix={}_LD --N={} --out={} ;".format(output,radius,output,samplesize,output)
    print cmd
    subprocess.call(cmd, shell = True)


def orig_summary(sample):
    gwas = {}
    with open(sample,'r') as f:
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
            chinese[rsid] = ':'.join([ref,alt,frequency])
    return chinese

def ldpred(output):
    with open('/share/data2/leon/gwas/ldpred/ldpred/' + output + '_LDpred-inf.txt','r') as f:
        inf = f.readlines()
    return inf



def final_result(inf,gwas,chinese,output):
    try:
        with open('/share/data2/leon/gwas/ldpred/{}/'.format(args.source) + output + '.adjusted.txt','w') as f:
            f.write('Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbJ\tbJ_se\tpJ\tLD_r\n')
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
                        f.write(str(1-float(chinese[rsid].split(':')[-1])) + '\t')
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

    except:
        cmd = "rm /share/data2/leon/gwas/ldpred/{}/{}.adjusted.txt".format(args.source,output)
        subprocess.call(cmd, shell = True)
        cmd = 'echo {} >>dirty_data'.format(output)
        print 'final result error'
        subprocess.call(cmd, shell = True)



if __name__ == '__main__':
    args = get_args()
    sample = args.input
    output = sample.split('/')[-1]
    source = args.source
    radius = args.radius
    if source == 'study':
        study_info = pd.read_table('/share/data2/leon/gwas/ldpred/study.csv',header = 0,sep = ',')
        samplesize = study_info[study_info['dataset'] == output]['n'].values[0]
    elif source == 'center':
        center_info = pd.read_table('/share/data2/leon/gwas/ldpred/center.csv',header = 0,sep = ',')
        samplesize = center_info[center_info['dataset'] == output]['Sample_size'].values[0]
    else:
        print 'only support share or gwas center source'
        sys.exit()
    #hsdf5(sample,output,samplesize)
    adjusted(output,samplesize,radius)
    gwas = orig_summary(sample)
    chinese = freq()
    inf = ldpred(output)
    final_result(inf,gwas,chinese,output)
    cmd = 'cp /share/data2/leon/gwas/ldpred/{}/{}.adjusted.txt /share/data2/leon/gwas/ldpred/step1/'.format(args.source,output)
    subprocess.call(cmd , shell = True)
