#!/usr/bin/python3
"""
optimized script for modeling the consequences of partial apomixis
each organism has mitochondrial DNA marker inherited from one of the parents 
and a set of microsats each evolving and segregating independently according simple 
SMM model by Kimura/Ohta 1973
==================================
	FULLY SEXUAL CASE
==================================
"""
import sys
import random
import time
import string
import numpy as np
import pandas as pd

version = "0.0.1"
start_time = time.time()
chars = ['A', 'C', 'G', 'T']
nmeChars = ['a', 'b', 'c', 't', 'p', 'l', 'u', 'x', '1', '2', '0', 'f', 'h', 'm']
ms_allel = [0b1, 0b10, 0b100, 0b1000]
homoGenotype = []
heteGenotype = []
for i in range(4):
	homoGenotype.append(ms_allel[i]+ms_allel[i])
	for j in range(i):
		heteGenotype.append(ms_allel[i]+ms_allel[j])
#-------------------------------------------------------------------------------------
lseq = 500
ms_num = 4 		#число локусов микросатов
ms_nallel = 4 	#число аллелей на локус
nGadov = 5000
nSteps = 10000
sampleSize = 20
nDNAMute = 500  #это число оснований, которые надо поменять в данном поколении
nMSMute = 10  #это число микросателлитных локусов, которые надо поменять в данном поколении
#-------------------------------------------------------------------------------------

def mkName():
	res = ""
	for i in range(4):
		ch = random.choice(nmeChars)
		res += ch
	return res

def sampleDNA(gg,nOTU,to_file=False):
	r = gg.sample(nOTU, axis=0, replace=False)
	seqs = []
	for rec in r:
		line = list(r["Seq"])
		l = ''
		#print(line)
		#for i in line:
			#print('>') 
		
	seqs = pd.DataFrame(r,columns=['Seq'])
	nmes = []
	sq = []
	count = 0
	for l in range(len(r)):
		nmes.append(str(count+1))
		count += 1
	res = pd.DataFrame([{'OTU':nmes,'Seq':seqs}])
	if(to_file):
		fout=open('sample.fas','w')
		count = 0
		for l in range(nOTU):
			q = ''
			for cc in r.iloc[l][0]:
				q += chars[cc]
			fout.write('>Seq_{}\n{}\n'.format(l,q))
			count += 1
		fout.close()
	return res
	
def sampleMSAT(gg,nOTU, to_file=False):
	r = gg.sample(nOTU, axis=0, replace=False)
	
def newGad(l,Nl,Na):
	#mtSeq = ''.join(np.random.choice(chars,l,replace=True))
	mtSeq = np.random.randint(0,4, l)
	#print(mtSeq)
	msats1 = []
	msats2 = []
	for i in range(Nl):
		set1 = ms_allel[np.random.randint(Na)]
		set2 = ms_allel[np.random.randint(Na)]
		msats1.append(set1)
		msats2.append(set2)
	gad = pd.DataFrame([{'Seq':mtSeq,'ms1':msats1,'ms2':msats2,'gender':0,'Lineage':0,'MSmute':0}])
	return gad

def cloneGad(g):
	return g.copy()

'''
JC model for mtDNA
'''
def mutate_DNA(c):
	y=c
	while(y==c):
		y=np.random.randint(0,4,1)
	z = y
	return(z[0])

#надо еще микросаты на 2 кучки разделить, отдает парой кассет
def mutate_ms(m1):
	res = np.random.choice(ms_allel)
	return(res)

def countUnique(g):
	#print(g["Lineage"])
	return(len(pd.unique(pd.Series(g["Lineage"]))))
	#return
	
def countMSatInfo(g):
	return(g.pivot_table(index=['Lineage'], aggfunc ='size'))
	
def main():
	fn = 'raw_'+mkName()+'.tab'
	inf = open(fn,'w')
	inf.write("\tstep\ta0\ta1\tb0\tb1\tc0\tc1\td0\td1\tnLineages\n")
	G0 = newGad(lseq, ms_num, ms_nallel)
	#G0 - основатель. The rest are random thus follow H-W rule
	Generation1=pd.DataFrame(columns=['Seq','ms1','ms2','gender','Lineage','MSmute'])	
	'''
	HW finding population
	'''
	for i in range(nGadov):
		Generation1 = pd.concat([Generation1,newGad(lseq, ms_num, ms_nallel)],axis=0,ignore_index=True)
	for i in range(nGadov):	
		Generation1.loc[i,"Lineage"] = i
		
	
	
	#Generation1 = Generation0.sample(nGadov,axis = 0,replace=True)
	MS_rg = nGadov*ms_num #ms array
	nLines = nGadov  #Число предковых линий исходно равно числу гадов
	for step in range(nSteps):
		Generation2 = Generation1.sample(nGadov, axis=0, replace=True,ignore_index=True)
		for i in range(nGadov):
			for ii in range(ms_num):
				tmp = Generation2.loc[i,"ms2"][ii]
				j=np.random.randint(nGadov)
				Generation2.loc[i,'ms2'][ii]=Generation2.loc[j,'ms2'][ii]
				Generation2.loc[j,'ms2'][ii]=tmp
		if ((step % 20)==0):
			print('Step # ',step)
			nLin=countUnique(Generation1)
			for j in range(nGadov):
			
				inf.write("{}\t{}\t{}\t{}\t".format(j,step,Generation2.loc[j,'ms1'][0],Generation2.loc[j,'ms2'][0]))
				inf.write("{}\t{}\t".format(Generation2.loc[j,'ms1'][1],Generation2.loc[j,'ms2'][1]))
				inf.write("{}\t{}\t".format(Generation2.loc[j,'ms1'][2],Generation2.loc[j,'ms2'][2]))
				inf.write("{}\t{}\t".format(Generation2.loc[j,'ms1'][3],Generation2.loc[j,'ms2'][3]))
				inf.write(str(nLin)+"\n")
			#print(step,countMSatInfo(Generation1))
		whoGoes = np.random.randint(0,MS_rg,nMSMute)
		nGad,nLoc = np.divmod(whoGoes, ms_num)
		
		for i in range(nMSMute):
			x = nGad[i]
			y = nLoc[i]
			#print(x,y)
			yy = np.random.randint(2)
			if(yy % 2 == 0):
				zz =Generation2.loc[x,'ms1'][y]
				Generation2.loc[x,'ms1'][y] = mutate_ms(zz)
				nLines += 1
				Generation2.loc[x,'Lineage'] = nLines
			else:
				zz =Generation2.loc[x,'ms2'][y]
				Generation2.loc[x,'ms2'][y] = mutate_ms(zz)
				nLines += 1
				Generation2.loc[x,'Lineage'] = nLines
	
			
#		D_mute = np.random.randint(0, rg, nDNAMute)
		#print(rg,D_mute)
#		a,b = np.divmod(D_mute, lseq)
#		for x in range(nDNAMute):
			# lt = Generation2.iloc[a[x]][0][b[x]]
# 			zt = mutate_DNA(lt)
# 			Generation2.iloc[a[x]][0][b[x]]=zt
			#print(lt,'->',Generation2.iloc[a[x]][0][b[x]])
		Generation1 = Generation2.copy()
		#print(np.sum(Generation1.index))
		
	#sampleDNA(Generation2, sampleSize, to_file=True)
	#Generation1.to_csv('Generation.csv')

	print('hoh')
	inf.close()

if __name__ == "__main__":
	main()
	print("--- %s seconds ---" % (time.time() - start_time))
