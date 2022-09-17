import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

os.getcwd()
os.chdir('C:/Users/nitis/Desktop/Davis RNAseq DEG')

df = pd.read_csv('dataset.csv', header=0)
df = df.set_index('GeneSymbol')

header_list = ["Name"]
#list = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list1 = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list2 = pd.read_csv('Heatmap_for_R21/gene_sets.gmt', names=header_list)
#rownames = list['Name'].to_list()
rownames = list1.append(list2) 
rownames = rownames['Name'].to_list()
my_list = df.columns.values.tolist()
colnames = my_list[2:25]
df = df[colnames]
df = df.loc[rownames]
df = df.T
############### Correlation analysis ###################
corr_df =  df.corr(method='pearson')
#corr_df = corr_df.dropna()
df_lt = corr_df.where(np.tril(np.ones(corr_df.shape)).astype(np.bool))
df_lt = df_lt.dropna(axis = 0, how = 'all')
df_lt = df_lt.dropna(axis = 1, how = 'all')

sns.set(style="white", font_scale=1)
f, ax = plt.subplots(figsize=(20, 20))
hmap=sns.heatmap(df_lt,cmap='coolwarm', vmin=-1, vmax=1, linewidths=0.5, annot=True, annot_kws={"size":5})
# display first few rows/columns of correlation matrix using iloc fucntion in Pandas
#hmap=sns.heatmap(df_lt,cmap="Spectral", annot=True)
hmap.figure.savefig("Correlation_Heatmap_Lower_Triangle.pdf",
                    format='pdf',
                    dpi=500)

#################################################
### Hippo signaling in high responder group #####
df = pd.read_csv('dataset_OLD_Highresp.csv', header=0)
df = df.rename(columns={'Unnamed: 0': 'Gene Symbol'})
df = df.set_index('Gene Symbol')

header_list = ["Name"]
#list = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list1 = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list2 = pd.read_csv('Heatmap_for_R21/gene_sets.gmt', names=header_list)
#rownames = list['Name'].to_list()
rownames = list1.append(list2) 
rownames = rownames['Name'].to_list()
my_list = df.columns.values.tolist()
colnames = my_list[1:7]
df = df[colnames]
df = df.loc[rownames]
df = df.T
############### Correlation analysis ###################
corr_df =  df.corr(method='pearson')
#corr_df = corr_df.dropna()
df_lt = corr_df.where(np.tril(np.ones(corr_df.shape)).astype(np.bool))
df_lt = df_lt.dropna(axis = 0, how = 'all')
df_lt = df_lt.dropna(axis = 1, how = 'all')

sns.set(style="white", font_scale=1)
f, ax = plt.subplots(figsize=(20, 20))
hmap=sns.heatmap(df_lt,cmap='coolwarm', vmin=-1, vmax=1, linewidths=0.5, annot=True, annot_kws={"size":5})
# display first few rows/columns of correlation matrix using iloc fucntion in Pandas
#hmap=sns.heatmap(df_lt,cmap="Spectral", annot=True)
hmap.figure.savefig("Correlation_Hippo_Old_HighResp.pdf",
                    format='pdf',
                    dpi=500)

#########################################################
##### Hippo signaling in high non-responder group #######
df = pd.read_csv('dataset_OLD_Nonresp.csv', header=0)
df = df.rename(columns={'Unnamed: 0': 'Gene Symbol'})
df = df.set_index('Gene Symbol')

header_list = ["Name"]
#list = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list1 = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list2 = pd.read_csv('Heatmap_for_R21/gene_sets.gmt', names=header_list)
#rownames = list['Name'].to_list()
rownames = list1.append(list2) 
rownames = rownames['Name'].to_list()
my_list = df.columns.values.tolist()
colnames = my_list[1:7]
df = df[colnames]
df = df.loc[rownames]
df = df.T
############### Correlation analysis ###################
corr_df =  df.corr(method='pearson')
#corr_df = corr_df.dropna()
df_lt = corr_df.where(np.tril(np.ones(corr_df.shape)).astype(np.bool))
df_lt = df_lt.dropna(axis = 0, how = 'all')
df_lt = df_lt.dropna(axis = 1, how = 'all')

sns.set(style="white", font_scale=1)
f, ax = plt.subplots(figsize=(20, 20))
hmap=sns.heatmap(df_lt,cmap='coolwarm', vmin=-1, vmax=1, linewidths=0.5, annot=True, annot_kws={"size":5})
# display first few rows/columns of correlation matrix using iloc fucntion in Pandas
#hmap=sns.heatmap(df_lt,cmap="Spectral", annot=True)
hmap.figure.savefig("Correlation_Hippo_NonResp.pdf",
                    format='pdf',
                    dpi=500)


#########################################################
##### Hippo signaling in high responder & non-responder group #######
df = pd.read_csv('dataset_Old_Resp_NonResp.csv', header=0)
df = df.rename(columns={'Unnamed: 0': 'Gene Symbol'})
df = df.set_index('Gene Symbol')

header_list = ["Name"]
#list = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list1 = pd.read_csv('CHOLESTEROL BIOSYNTHESIS REACTOME.txt', names=header_list)
list2 = pd.read_csv('Heatmap_for_R21/gene_sets.gmt', names=header_list)
#rownames = list['Name'].to_list()
rownames = list1.append(list2) 
rownames = rownames['Name'].to_list()
my_list = df.columns.values.tolist()
colnames = my_list[1:7]
df = df[colnames]
df = df.loc[rownames]
df = df.T
############### Correlation analysis ###################
corr_df =  df.corr(method='pearson')
#corr_df = corr_df.dropna()
df_lt = corr_df.where(np.tril(np.ones(corr_df.shape)).astype(np.bool))
df_lt = df_lt.dropna(axis = 0, how = 'all')
df_lt = df_lt.dropna(axis = 1, how = 'all')

sns.set(style="white", font_scale=1)
f, ax = plt.subplots(figsize=(20, 20))
hmap=sns.heatmap(df_lt,cmap='coolwarm', vmin=-1, vmax=1, linewidths=0.5, annot=True, annot_kws={"size":5})
# display first few rows/columns of correlation matrix using iloc fucntion in Pandas
#hmap=sns.heatmap(df_lt,cmap="Spectral", annot=True)
hmap.figure.savefig("Correlation_Hippo_Resp_and_NonResp.pdf",
                    format='pdf',
                    dpi=500)

