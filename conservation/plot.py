import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import stats
import math
def plot_av(files):
    hd=[]
    ho=[]
    f=csv.reader(open(files, 'rb'))
    for row in f:
       x=row[0].split()
       hd.append(x[5])
       ho.append(x[6])
    del hd[0]
    del ho[0]
    
    rem_ent=[]
    
    for n in range(len(ho)):
	if ho[n]=='*':
		
		rem_ent.append(n)
    for n in range(len(hd)-1):
	
	if hd[n]=='*':
		rem_ent.append(n)
    
    offset=0   
    for n in sorted(rem_ent):
	del hd[n-offset]
	del ho[n-offset]
	offset+=1
  
    
    hd=np.array(hd,dtype=np.float32)
    ho=np.array(ho,dtype=np.float32)
    s,p=stats.shapiro(hd)
    
    t,p=stats.ttest_rel(hd,ho)
    
    print t,p
    print s,p
    plt.plot(hd,ho, 'g.')
    plt.xlabel('<Hd>')
    plt.ylabel('<Ho>')
    #plt.text(3.0,0.5, 't='+str(round(t,2)))
    #plt.text(3.0,0.3, 'p-value='+str(round(p,2)))
    plt.savefig('Hd_vs_Ho.pdf')
    plt.show()
    



def plot(files):
    hom=[]
    ent=[]
    f=csv.reader(open(files,'rb'))
    for row in f:
        x=row[0].split()
    	hom.append(x[4])
	ent.append(x[7])
    del hom[0]
    del ent[0]
    
    #t-test#
    
    hom=np.array(hom,dtype=np.float32)
    ent=np.array(ent,dtype=np.float32)
   
    two_sample=stats.ttest_rel(hom,ent)
    print two_sample
    plt.plot(hom,ent, 'r.')
    plt.xlabel('Number homologs')
    plt.ylabel('<Entropy>=<Hd>/<Ho>')
    plt.savefig('nhom_vs_entropy.pdf')
    plt.show()
   
def lenvsreg(files):
	dis_len=[]
	dis_ent=[]
	ord_len=[]
	ord_ent=[]

	f=csv.reader(open(files,'rb'))
	for row in f:
		x=row[0].split()
		dis_len.append(x[1])
		dis_ent.append(x[2])
		ord_len.append(x[3])
		ord_ent.append(x[4])
	del dis_len[0]
	del dis_ent[0]
	del ord_len[0]
	del ord_ent[0]
	
	dis_ent=[x for x in dis_ent if x!= 'None']
	dis_len=[x for x in dis_len if x!= 'None']
	ord_ent=[x for x in ord_ent if x!= 'None']
	ord_len=[x for x in ord_len if x!= 'None']
	
	plt.plot(dis_len,dis_ent, 'r.')
	plt.plot(ord_len,ord_ent, 'g.')
	plt.show()
	
def histogram(files):
	ent=[]
	f=csv.reader(open(files, 'rb'))
		
	for row in f:
        	x=row[0].split()
		ent.append(x[7])
	del ent[0]
	
	ent=np.array(ent,dtype=np.float32)
	
	plt.hist(ent, bins=30)
	plt.ylabel('Frequencies')
    	plt.xlabel('<Entropy>=<Hd>/<Ho>')
	
	plt.savefig('Histogram.pdf')
	plt.show()
def ttest(files):
	ttest=[]
	f=csv.reader(open(files, 'rb'))
		
	for row in f:
        	x=row[0].split()
		ttest.append(x[3])
	del ttest[0]
	rem_tt=[]
	for n in range(len(ttest)):
		if ttest[n]=='*':
			rem_tt.append(n)
    
    	offset=0   
   	for n in sorted(rem_tt):
		del ttest[n-offset]
		offset+=1
	ttest=np.array(ttest,dtype=np.float32)
	plt.title('t-test')
	plt.hist(ttest, bins=30)
	plt.savefig('ttest.pdf')
	plt.show()
	
ttest('/home/kettina/Scrivania/ttest.txt')
plot_av('/home/kettina/Scrivania/result.txt')
plot('/home/kettina/Scrivania/result.txt')
lenvsreg('/home/kettina/Scrivania/res_len.txt')
histogram('/home/kettina/Scrivania/result.txt')
