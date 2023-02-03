'''Script by Bethany Moore
PURPOSE: calculate percent overlap of clusters 
INPUT:
	REQUIRED:
	-bin binary matrix with clusters you are interested in
	-mr mr output with all clusters
	OPTIONAL
OUTPUT:
	dataframe of your cluster overlap percent (_clusteroverlap.txt)
	percentiles of all cluster overlap (_allcluster_percentiles.txt)

'''

import sys, os
import numpy as np
import pandas as pd

for i in range (1,len(sys.argv),2):
	if sys.argv[i] == '-bin': 
		bin_file = sys.argv[i+1]
	if sys.argv[i] == '-mr': 
		mr_file = sys.argv[i+1]

def get_bin(bin_file,D):
	inp= open(bin_file,'r')
	header= inp.readline()
	clusterlist= header.strip().split('\t')[1:]
	print(clusterlist)
	for line in inp:
		L= line.strip().split('\t')
		gene= L[0]
		data= L[1:]
		#print(data)
		for i in range(len(data)):
			j= data[i]
			if int(j) == 1:
				clust= clusterlist[i]
				if clust not in D:
					D[clust]=[gene]
				else:
					D[clust].append(gene)
			else:
				pass
	return D

def get_overlap(D,D2):
	key_list= list(D.keys())
	#print(key_list)
	for key in key_list:
		genelist1=D[key]
		for key2 in key_list:
			if key2 != key:
				genelist2=D[key2]
				set1 = set(genelist1) #get unique set for each gene list
				set2 = set(genelist2)
				overlap = list(set1 & set2) #get overlap
				universe = list(set1 | set2) # get union (all from both sets)
				result1 = float(len(overlap)) / len(set1) * 100 #calculate percent
				result2 = float(len(overlap)) / len(set2) * 100
				result3 = float(len(overlap)) / len(universe) * 100
				if (key2,key) not in D2.keys():
					D2[(key,key2)]=[result1,result2,result3]
		#key_list.remove(key)
		#print(key_list)
	return D2

def get_all_clusters(mr_file,D): # get all clusters in data- this is output from mr scripts
	inp= open(mr_file,'r')
	for line in inp:
		L= line.strip().split('\t')
		clust= L[0]
		data= L[3]
		genes= data.split(' ')
		D[clust]=genes
	return D

D={}
D2={}
D3={}
D4={}
# get clusters that you want from binary matrix
clust1_D= get_bin(bin_file, D)
# get overlap for these clusters
overlap1_D= get_overlap(clust1_D,D2)
# write out overlap for these clusters
df_overlap1=pd.DataFrame.from_dict(overlap1_D, orient='index') # convert to dataframe
print(df_overlap1)
df_overlap1 = df_overlap1.sort_index(ascending=True)
print(df_overlap1)
df_overlap1.to_csv(path_or_buf=str(bin_file)+"_clusteroverlap.txt", sep='\t',header=True, index=True)
# get all clusters in dataset
clust2_D= get_all_clusters(mr_file,D3)
# calculate overlap for all clusters
overlap2_D= get_overlap(clust2_D,D4)
# calculate percentiles for all cluster overlap
df_overlap=pd.DataFrame.from_dict(overlap2_D) # convert to dataframe
df1= df_overlap.quantile(q=0.95, axis=1) # get percentile
df2= df_overlap.quantile(q=0.99, axis=1)
df3= df_overlap.quantile(q=0.995, axis=1)
# combine dataframes, write out
frames = [df1, df2, df3]
result = pd.concat(frames, axis=1)
print(result)
result.to_csv(path_or_buf=str(bin_file)+"_allcluster_percentiles.txt", sep='\t',header=True, index=True)

# 				