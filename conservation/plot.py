
import matplotlib.pyplot as plt
import csv

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
  
    print hd
    print ho
    plt.plot(hd,ho, 'go')
    plt.xlabel('<Hd>')
    plt.ylabel('<Ho>')
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
    print hom
    print ent
    plt.plot(hom,ent, 'ro')
    plt.xlabel('Number homologs')
    plt.ylabel('<Entropy>')
    plt.show()
'''
def plot_species(files,spec):
    species={}
    entropies={}
    f=open(spec,'r')
    f=f.readlines()
    f=[i.split(',') for i in f]
    for line in f:

        species[line[0]]=line[3]
    
    f=csv.reader(open(files,'rb'))
    for row in f:
        x=row[0].split()
    	entropies[x[0]]=x[7]
 
    
        
    del species['disprotID']
    spec=[]
    for x in species:
	spec.append(species[x])
    spec=set(spec)
    
	#if species[x]=='Homo sapiens':
	#	hom_sap.appendentropies[x]
plot_species('/home/enrichetta/Desktop/prova.txt', '/media/data/dataset_oxana/proteins.txt')
'''
plot_av('/home/enrichetta/Desktop/prova.txt')
plot('/home/enrichetta/Desktop/prova.txt')
