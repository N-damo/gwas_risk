#!/share/data3/lianlin/soft/python/bin/python2.7
import pandas as pd 
import os

def merge_file():
    df = pd.read_table('/share/data2/leon/report/info.txt', header = 0, sep = '\t',encoding='utf_8')
    first = df.loc[0,'id'] 
    assert os.path.exists(first + '.csv'),'error {} does not exist, please check it'.format(first + '.csv')
    merge = pd.read_table(first + '.csv', header = 0, sep = ',')
    print len(merge)
    for i in df['id'][1:].values:
        merge2 = pd.read_table(i + '.csv', header = 0, sep = ',')
        merge = pd.merge(merge, merge2, on = 'sample')
    merge=merge.set_index('sample')
    merge.loc['avg']=merge.mean()
    merge.loc['std']=merge.std()
    merge=merge.reset_index()
    merge.rename(columns={'sample':'index'},inplace=True)
    # for i in range(len(merge)):
    #     merge.loc[i,]
    return merge

if __name__ == '__main__':
    merge = merge_file()
    merge.to_csv('polygenic_risk.csv', header = True, index = False, encoding='utf_8')

