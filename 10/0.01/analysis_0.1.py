#!/usr/bin/python3

import sys
import time
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

dta_file1 ='raw_1ptm.tab'
start_time = time.time()

def count_zeros(df):
	y = len(df.columns)
	x = len(df.index)
	res = np.zeros(x)
	for i in range(x):
		line = df.iloc[i]
		res[i] = np.count_nonzero(line == 0)
	return(res)

def main():
	
	dta = pd.read_csv(dta_file1,sep='\t')
	ngadov = int(len(dta.index)/len(dta['step'].unique()))
	print(ngadov)
	allels = [1,2,4,8]
	doubles = []
	e_bins = []
	for i in range(4):
		for j in range(4):
			doubles.append(16*allels[i]+allels[j])
			e_bins.append(16*allels[i]+allels[j])
	e_bins.append(200)
	
	a0b0 = 16*dta["a0"]+dta["b0"]
	a1b1 = 16*dta["a1"]+dta["b1"]
	a0c0 = 16*dta["a0"]+dta["c0"]
	a1c1 = 16*dta["a1"]+dta["c1"]
	a0d0 = 16*dta["a0"]+dta["d0"]
	a1d1 = 16*dta["a1"]+dta["d1"]
	b0c0 = 16*dta["b0"]+dta["c0"]
	b1c1 = 16*dta["b1"]+dta["c1"]
	b0d0 = 16*dta["b0"]+dta["d0"]
	b1d1 = 16*dta["b1"]+dta["d1"]
	c0d0 = 16*dta["c0"]+dta["d0"]
	c1d1 = 16*dta["c1"]+dta["d1"]
	stps = dta['step'].unique()
	stpss = dta['step']
	ndata = len(dta['step'])
	nlines = ndata*16
	haplo_names = ['a0b0','a1b1','a0c0','a1c1','a0d0','a1d1','b0c0','b1c1','b0d0','b1d1','c0d0','c1d1']
	haplos = pd.DataFrame({'N':np.arange(ndata),'step':dta["step"],'a0b0':a0b0, 'a1b1':a1b1,'a0c0':a0c0,
	'a1c1':a1c1,'a0d0':a0d0,'a1d1':a1d1,'b0c0':b0c0,'b1c1':b1c1,'b0d0':b0d0,'b1d1':b1d1,'c0d0':c0d0,'c1d1':c1d1})
	
	print('here')
	print('doing density matrices')
	d2_a0b0 = []
	d2_a1b1 = []
	d2_a0c0 = []
	d2_a1c1 = []
	d2_a0d0 = []	
	d2_a1d1 = []
	d2_b0c0 = []
	d2_b1c1 = []
	d2_b0d0 = []			
	d2_b1d1 = []
	d2_c0d0 = []
	d2_c1d1 = []
	metrics=pd.DataFrame(columns=['average'])
	for shag in stps:
		generation = haplos.set_index("step").loc[shag]
		d2_a0b0.append(list(np.histogram(generation['a0b0'],bins=e_bins)[0]))
		d2_a1b1.append(list(np.histogram(generation['a1b1'],bins=e_bins)[0]))
		d2_a0c0.append(list(np.histogram(generation['a0c0'],bins=e_bins)[0]))
		d2_a1c1.append(list(np.histogram(generation['a1c1'],bins=e_bins)[0]))
		d2_a0d0.append(list(np.histogram(generation['a0d0'],bins=e_bins)[0]))
		d2_a1d1.append(list(np.histogram(generation['a1d1'],bins=e_bins)[0]))
		d2_b0c0.append(list(np.histogram(generation['b0c0'],bins=e_bins)[0]))
		d2_b1c1.append(list(np.histogram(generation['b1c1'],bins=e_bins)[0]))
		d2_b0d0.append(list(np.histogram(generation['b0d0'],bins=e_bins)[0]))
		d2_b1d1.append(list(np.histogram(generation['b1d1'],bins=e_bins)[0]))
		d2_c0d0.append(list(np.histogram(generation['c0d0'],bins=e_bins)[0]))
		d2_c1d1.append(list(np.histogram(generation['c1d1'],bins=e_bins)[0]))
#counting zeros in a line		
	dst0 = list(count_zeros(pd.DataFrame(d2_a0b0,columns = doubles,index=stps)))
	dst1 = list(count_zeros(pd.DataFrame(d2_a1b1,columns = doubles,index=stps)))
	dst2 = list(count_zeros(pd.DataFrame(d2_a0c0,columns = doubles,index=stps)))
	dst3 = list(count_zeros(pd.DataFrame(d2_a1c1,columns = doubles,index=stps)))
	dst4 = list(count_zeros(pd.DataFrame(d2_a0d0,columns = doubles,index=stps)))
	dst5 = list(count_zeros(pd.DataFrame(d2_a1d1,columns = doubles,index=stps)))
	dst6 = list(count_zeros(pd.DataFrame(d2_b0c0,columns = doubles,index=stps)))
	dst7 = list(count_zeros(pd.DataFrame(d2_b1c1,columns = doubles,index=stps)))
	dst8 = list(count_zeros(pd.DataFrame(d2_b0d0,columns = doubles,index=stps)))
	dst9 = list(count_zeros(pd.DataFrame(d2_b1d1,columns = doubles,index=stps)))
	dst10 = list(count_zeros(pd.DataFrame(d2_c0d0,columns = doubles,index=stps)))
	dst11 = list(count_zeros(pd.DataFrame(d2_c1d1,columns = doubles,index=stps)))
	
#make a csv file for R	
	count_0=dst0+dst1+dst2+dst3+dst4+dst5+dst6+dst7+dst8+dst9+dst10+dst11
	generation = list(np.tile(np.arange(len(stps)),12))
	haplotype = list(np.repeat(haplo_names,len(stps)))
	#print(len(count_0),len(generation),len(haplotype))
	zrs = pd.DataFrame({'generation':generation,'haplotype':haplotype,'count_0':count_0})
	zrs.to_csv('zeros0.01.tab')
	
	sns.heatmap(dst0)
	plt.show()

if __name__ == "__main__":
	main()
	print("--- %s seconds ---" % (time.time() - start_time))