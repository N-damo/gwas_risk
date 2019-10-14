#!/usr/bin/env python3
# coding:utf-8

import pandas as pd
import numpy as np
import morbidity
from scipy.stats import norm
import math
import json
from collections import defaultdict
import sys
#import pysnooper



class Disease(morbidity.Beta):
    #trait_list='/share/data3/lianlin/soft/bin/gwas/risk/data/test.csv'
    def __init__(self, vcf, all_summary_data, beta_avg_std,trait_list):
        super(Disease, self).__init__(vcf, all_summary_data)
        self.trait_list=trait_list
        self.cnmorbidity_table=self.disired_table()
        self.beta_avg_std = beta_avg_std


    def disired_table(self):
        #挑选感兴趣的疾病，这些疾病具有人群患病率
        cnmorbidity_table = pd.read_csv(self.trait_list, sep=',', header=0, encoding='gbk')
        need_columns=['sex','morbidity1','morbidity2','morbidity3','morbidity4','gene_factor','English','Chinese','basic_data_name']
        #available_columns=cnmorbidity_table.columns()
        try:
            cnmorbidity_table=cnmorbidity_table[need_columns]
        except KeyError:
            print('we need these columns,plese make sure\n{}'.format(need_columns))
 
        trait=list(cnmorbidity_table['basic_data_name'])
        self.all_dataset={key:value for key,value in self.all_dataset.items() if key in trait}
        return cnmorbidity_table



    def sub_cnmorbidity(self, gender, dataset_name):
        dataset = self.cnmorbidity_table[self.cnmorbidity_table['basic_data_name'] == dataset_name]
        phenotype=dataset['English']
        common_pop_list={}
        for phe in phenotype:
            dataset_sub=dataset[dataset['English'] == phe]
            Chinese_phe=dataset_sub['Chinese'].values[0]
            #print(dataset_name)
            trait_include_male_and_female = int(dataset_sub['sex'])
            gene_factor=dataset_sub['gene_factor']
            # male=morbidity1,morbidity3
            # female=morbidity2,morbidity4
            # sex=0,only male
            # sex=1,only female
            #sex=2,male and female4
            if trait_include_male_and_female == 2:
                if gender == 0:
                    common_pop1 = dataset_sub['morbidity1'].values[0]
                    common_pop2 = dataset_sub['morbidity3'].values[0]
                else:
                    common_pop1 = dataset_sub['morbidity2'].values[0]
                    common_pop2 = dataset_sub['morbidity4'].values[0]

            elif trait_include_male_and_female == 0:
                if gender == 0:
                    common_pop1 = dataset_sub['morbidity1'].values[0]
                    common_pop2 = dataset_sub['morbidity3'].values[0]
                else:
                    common_pop1 = 0
                    common_pop2 = 0

            elif trait_include_male_and_female == 1:
                if gender == 1:
                    common_pop1 = dataset_sub['morbidity2'].values[0]
                    common_pop2 = dataset_sub['morbidity4'].values[0]
                else:
                    common_pop1 = 0
                    common_pop2 = 0

            if np.isnan(common_pop1):
                common_pop = common_pop2
            else:
                common_pop = common_pop1
            common_pop_list[phe]=[float(common_pop),float(gene_factor),Chinese_phe]
        return common_pop_list

    def load_avg_std(self):
        with open(self.beta_avg_std, 'r') as f:
            for i in f:
                line = i.strip()
                beta_avg_std = json.loads(line)

        return beta_avg_std
        
    #@pysnooper.snoop('./test.log')
    def risk_level(self):
        report = defaultdict(dict)
        sample_beta = self.parse_variant()
        sample_gender = self.infer_gender()
        beta_avg_std = self.load_avg_std()
        for dataset in sample_beta:
            for sample in sample_beta[dataset]:
                beta = sample_beta[dataset][sample]
                gender = sample_gender[sample]
                common_pop_list = self.sub_cnmorbidity(gender, dataset)
                for phe in common_pop_list:
                    cnmorbidity,gene_factor,Chinese_phe=common_pop_list[phe]
                    if cnmorbidity == 0:
                        continue
                    else:
                        avg, std = beta_avg_std[dataset]
                        t = norm.ppf(1-cnmorbidity)  # 人群发病率对应beta值
                        new_beta = (beta - avg)/std
                        x = new_beta-1
                        _ = 1-norm.cdf(t-x)
                        individual_morbidity = _*gene_factor  # average h = 20%
                        times = individual_morbidity/cnmorbidity
                        if times > 1:
                            original_risk_value = 1+math.log(times, math.e**2)
                        else:
                            original_risk_value = times
                        # if times >1,risk_value may still lower than 1
                        # (1+math.log(math.e, math.e**2))=1.5
                        risk_value = original_risk_value / \
                            (1+math.log(math.e, math.e**2))
                        suggest = self.report_value(
                            risk_value, cnmorbidity, individual_morbidity)
                        #Chinese_phe=self.cnmorbidity_table[self.cnmorbidity_table['English'] == phe]['Chinese'].values[0]
                        report[sample][phe+'('+Chinese_phe+')']=[gender,cnmorbidity,individual_morbidity, risk_value, suggest,dataset]
                        #print(sample,gender, phe,Chinese_phe,cnmorbidity,individual_morbidity, risk_value, suggest, sep='\t')
        return report

    def morbidity_level(self, morbidity):
        if morbidity >= 1e-1:
            level = 3
        elif 1e-3 <= morbidity < 1e-1:
            level = 2
        else:
            level = 1
        return level

    def individual_risk_value(self, risk_value):
        # (1+math.log(50*math.e, math.e**2))/(1+math.log(math.e, math.e**2)) ≈ 2.3
        if risk_value < 1:
            level = 1
        elif 1 <= risk_value < (1+math.log(50*math.e, math.e**2))/(1+math.log(math.e, math.e**2)):
            level = 2
        else:
            level = 3
        return level

    def report_value(self, risk_value, cnmorbidity, individual_morbidity):
        cal_cnmorbidity_level = self.morbidity_level(cnmorbidity)
        cal_individual_morbidity_level = self.morbidity_level(
            individual_morbidity)
        if risk_value == 1:
            if cal_individual_morbidity_level > cal_cnmorbidity_level:
                cal_individual_morbidity_level = cal_cnmorbidity_level
        individual_risk_value_level = self.individual_risk_value(risk_value)
        if cal_individual_morbidity_level == 3 and individual_risk_value_level == 3:
            return '重点关注'
        elif cal_individual_morbidity_level == 3 or individual_risk_value_level == 3:
            return '适度关注'
        elif individual_risk_value_level == 2 or cal_individual_morbidity_level == 2:
            return '略微关注'
        else:
            return '正常对待'

def json_out(report):
    with open('report_disease.json','w') as f:
        json.dump(report,f,ensure_ascii=False)



if __name__ == '__main__':
    vcf=sys.argv[1]
    # vcf = '/share/data4/deposit/VCF/gmb21.vcf.gz'
    trait_list='/share/data3/lianlin/soft/bin/gwas/risk/data/test.csv'
    all_summary_data = '/share/data3/lianlin/soft/bin/gwas/risk/data/trait_locus.json'
    beta_avg_std = '/share/data3/lianlin/soft/bin/gwas/risk/data/beta_avg_std.json'
    d = Disease(vcf, all_summary_data, beta_avg_std,trait_list)
    report=d.risk_level()
    df=pd.DataFrame.from_dict(report,orient='index').T
    df=df.reset_index()
    df.to_csv('report.csv',header=True,index=False,encoding='gbk')
    json_out(report)
