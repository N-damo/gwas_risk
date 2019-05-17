#!/usr/bin/python2.7
import gzip
import sys
import math
import re
import argparse
import subprocess
from argparse import ArgumentParser
def get_args():
    parser=ArgumentParser(description='gwas summary data format')
    parser.add_argument('-i','--input',help='input your file')
    args=parser.parse_args()
    return args

# def hg38():
#     with open('/share/data2/leon/chinese','r') as f:
#         l = f.realines()
#     return l

def is_number(num):
    pattern = re.compile(r'^[-+]?[-0-9]\d*\.\d*|[-+]?\.?[0-9]\d*$')
    result = pattern.match(num)
    if result:
        return True
    else:
        return False

def study(output,l):
    print output
    with open('/share/data2/leon/gwas/study/basic1/' + output,'w') as f:
        for i in l[1:]:
            h=i.strip().split('\t')
            beta = h[4]
            a1 = h[2]
            a2 = h[3]
            pvalue = float(h[6])
            rsid = h[0]
            if beta == '-' :
                print 'dirty data'
                cmd1 = 'echo {} >>dirty_data'.format(output)
                cmd2 = 'rm {} '.format(output)
                subprocess.call(cmd1, shell=True)
                subprocess.call(cmd2, shell=True)
                sys.exit()
            elif is_number(beta):
                if 0<pvalue<1:            
                    f.write(h[1].split(':')[0] + '\t')
                    f.write(rsid + '\t')
                    f.write(a1 + '\t')
                    f.write(a2 + '\t')
                    f.write(h[1].split(':')[1] +'\t')
                    f.write(str(math.e**float(beta)) + '\t')
                    f.write(str(pvalue) + '\n')
                    
                    



if __name__ == '__main__':
    args = get_args()
    sample = args.input
    with gzip.open(sample,'r') as f:
        l=f.readlines()
    sample = sample.split('/')[-1]
    output = sample.replace('.gz', '')
    # hg = hg38()
    study(output,l)