import os
import sys


files=os.listdir(sys.argv[1])
js_file=[f for f in files if f != 'cut_trait.json']

sge='#$ -N genotype\n#$ -pe smp 3\n#$ -q all.q\n#$ -cwd\nset -e\ncd /share/data2/leon/polygenic/beta\nsource ~/.bash_profile'


for js in js_file:
    with open(os.path.basename(js),'w') as f:
        f.write(sge+'\n')
        for i in range(8,22):
            cmd="python3 /share/data3/lianlin/soft/bin/gwas/risk/collect_beta.py /share/data4/deposit/VCF/gmb{}.vcf.gz /share/data2/leon/polygenic/sub/{} ;".format(i,js)
            f.write(cmd+'\n')









