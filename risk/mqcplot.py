#!/usr/bin/env python3
#coding:utf-8

import json
import os
import sys
from collections import defaultdict

"""
{
    "id": "custom_data_json_table",
    "section_name": "Custom JSON Table",
    "description": "This table is a self-contained JSON file. Hopefully the column order will be honoured.",
    "plot_type": "table",
    "pconfig": {
        "id": "custom_data_json_table_table",
        "title": "JSON table with ordered columns",
        "min": 0
    },
    "data": {
        "sample_1": {
            "1": 12,
            "2": 14,
            "3": 10,
            "4": 7,
            "5": 16
        },
        "sample_2": {
            "1": 9,
            "2": 11,
            "3": 15,
            "4": 18,
            "5": 21
        }
    }
}
"""


def nesteddict():
    return defaultdict(nesteddict)

    
class Multiqc(object):
    """ report[sample][phe+'('+Chinese_phe+')']=[gender,cnmorbidity,individual_morbidity, risk_value, suggest,dataset] """
    def __init__(self,disease:dict):
        self.samples=disease
        self.sample_gender = nesteddict()
        self.phe_cnmorbidity = nesteddict()
        self.sample_morbidity = nesteddict()
        self.sample_risk_value = nesteddict()
        self.sample_suggest = nesteddict()
        self.health_distribution = nesteddict()
        #self.special_view = nesteddict()
        for sample in self.samples:
            self.health_distribution[sample]['重点关注']=0
            self.health_distribution[sample]['适度关注']=0
            self.health_distribution[sample]['略微关注']=0
            self.health_distribution[sample]['正常对待']=0

        for sample in self.samples:
            for phe in self.samples[sample]:
                gender,cnmorbidity,individual_morbidity, risk_value, suggest,_ = self.samples[sample][phe]
                self.sample_gender[sample]['gender']=['male' if gender == 0 else 'female'][0]
                self.phe_cnmorbidity[phe]['cnmorbidity']=format(cnmorbidity,'.2E')
                self.sample_morbidity[phe][sample]=format(individual_morbidity,'.2E')
                self.sample_risk_value[phe][sample]=risk_value
                self.sample_suggest[phe][sample]=suggest
                self.health_distribution[sample][suggest] += 1


        self.gender_plot()
        self.cnmorbidity_plot()
        self.morbidity_plot()
        self.risk_value_plot()
        self.suggest_plot()
        self.stat_percent_plot()



    def out_json(self,plot:dict,name:str):
        with open(name+'_mqc.json','wt') as f:
            json.dump(plot,f,ensure_ascii=False)


    def gender_plot(self):
        sample_gender_plot = nesteddict()
        sample_gender_plot['id']="sample_gender_table"
        sample_gender_plot["section_name"]="Sample Gender Table"
        sample_gender_plot["description"]="sample dengder predicted by x chromosome snp hom rate."
        sample_gender_plot["plot_type"]="table"
        sample_gender_plot["pconfig"]={"id": "sample_gender_table_table","title": "Gender Predict"}
        sample_gender_plot["data"]=self.sample_gender
        return self.out_json(sample_gender_plot,'gender_plot')


    def cnmorbidity_plot(self):
        phe_cnmorbidity_plot = nesteddict()
        phe_cnmorbidity_plot['id']="disease_cnmorbidity_table"
        phe_cnmorbidity_plot["section_name"]="Disease Cnmorbidity Table"
        phe_cnmorbidity_plot["description"]="Some disease common morbidity collected from web source or self-predicted."
        phe_cnmorbidity_plot["plot_type"]="table"
        phe_cnmorbidity_plot["pconfig"]={"id": "disease_cnmorbidity_table_table","title": "disease common morbidity", "min": 0}
        phe_cnmorbidity_plot["data"]=self.phe_cnmorbidity
        return self.out_json(phe_cnmorbidity_plot,'cnmorbidity_plot')

    
    def morbidity_plot(self):
        sample_morbidity_plot = nesteddict()
        sample_morbidity_plot['id']="sample_morbidity_table"
        sample_morbidity_plot["section_name"]="Individual's Morbidity Table"
        sample_morbidity_plot["description"]="Individual's disease morbidity adjusted by normal distribution."
        sample_morbidity_plot["plot_type"]="table"
        sample_morbidity_plot["pconfig"]={"id": "sample_morbidity_table_table","title": "individual's morbidity", "min": 0}
        sample_morbidity_plot["data"]=self.sample_morbidity
        return self.out_json(sample_morbidity_plot,'sample_morbidity')
    
    def risk_value_plot(self):
        sample_risk_value_plot = nesteddict()
        sample_risk_value_plot['id']="sample_risk_value_table"
        sample_risk_value_plot["section_name"]="Risk Value Table"
        sample_risk_value_plot["description"]="This table is a self-contained JSON file. Hopefully the column order will be honoured."
        sample_risk_value_plot["plot_type"]="table"
        sample_risk_value_plot["pconfig"]={"id": "custom_data_json_table_table","title": "JSON table with ordered columns", "min": 0}
        sample_risk_value_plot["data"]=self.sample_risk_value
        return self.out_json(sample_risk_value_plot,'sample_risk')

    def suggest_plot(self):
        sample_suggest_plot = nesteddict()
        sample_suggest_plot['id']="sample_suggest_table"
        sample_suggest_plot["section_name"]="Individual's health suggest Table"
        sample_suggest_plot["description"]="some health suggest will be showed according to cnmorbidy,individual's morbidty and risk value."
        sample_suggest_plot["plot_type"]="table"
        sample_suggest_plot["pconfig"]={"id": "sample_suggest_table_table","title": "health suggest", "min": 0}
        sample_suggest_plot["data"]=self.sample_suggest  
        return self.out_json(sample_suggest_plot,'sample_suggest')   

    def stat_percent_plot(self):
        stat_percent_dis = nesteddict()
        stat_percent_dis['id']="stat_percent_dis"
        stat_percent_dis["section_name"]="Health Distribution"
        stat_percent_dis["description"]="Health suggestion distribution."
        stat_percent_dis["plot_type"]="bargraph"
        stat_percent_dis["pconfig"]={"id": "stat_percent_dis_dis","title": "health distribution", "min": 0}
        stat_percent_dis["data"]=self.health_distribution  
        return self.out_json(stat_percent_dis,'stat_percent')   





def load_json(report:'json'):
    with open(report,'rt') as f:
        for i in f:
            line=i.strip()
            report=json.loads(line)
    return report


if __name__ == '__main__':
    report = sys.argv[1]
    report = load_json(report)
    Multiqc(report)