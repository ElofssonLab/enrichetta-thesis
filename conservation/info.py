import csv
import urllib2
import os
import glob
from  goatools import obo_parser 



def ids(files):
	f=csv.reader(open(files,'rb'))
	global dic_id
	dic_id={}
	for row in f:
		x=row[0].split()
		dic_id[x[0]]=x[1]
	del dic_id['Disprot_ID']
	
def info(files):
	global ids
	ids=[]	
	f=csv.reader(open(files,'rb'))
	for row in f:
		x=row[0].split()
		if x[2]!='*':
			ids.append(x[2])
	del ids[0]	
def download():
	global el
	for el in ids:
		if el in dic_id.keys():
			print dic_id[el], el
			file=dic_id[el]+'.txt'
			url='http://www.uniprot.org/uniprot/'+dic_id[el]+'.txt'
			print url
			response=urllib2.urlopen(url)
			with open(os.path.join('/home/enrichetta/Project/download',file), 'w') as f:
				f.write(response.read())
def parse(files):
	global name
	global GOs,Pfam
	global sp,p
	global species
	global idss
	drive, path=os.path.splitdrive(files)
     	path, filename=os.path.split(path)
     	name=filename.split('.')[0]
	GOs=[]
	Pfam=[]
	species={}
	idss={}
	with open(files,'r') as f:
		f=f.readlines()
		for line in f:
			if 'GO' in line:
				GOs.append(line.split(';')[1].strip())
			if 'Pfam' in line:
				x=line.split(';')[1].strip()
				Pfam.append(x)
				print Pfam
			
					
			if line.startswith('OS'):
				oS= ' '.join(line.split()[1:])
				
				species[name]=oS
				
			if line.startswith('ID'):
				i=line.split()[1]
				idss[name]=i
			
				
				
				
	sp={'cellular_component':'C', 'biological_process':'P', 'molecular_function':'F'}
	
	p=obo_parser.GODag('/home/enrichetta/goatools/data/gene_ontology.1_2.obo')
def printparse():
	for key,value in dic_id.items():
		if value==name:
			print key,value, idss[name],Pfam[name],species[name]
def go():				
	
	
	for key,value in dic_id.items():
		if value==name:
			for l in GOs:
		
				print '{0:12} {1:12} {2:15} {3:5}{4:12}'.format(key,value,l, sp[p[l].namespace], p[l].name)
	


		

		



				
ids('/home/enrichetta/Project/enrichetta-thesis/dataset.txt')
info('/home/enrichetta/Project/enrichetta-thesis/conservation/cons2.txt')
#download()
path='/home/enrichetta/Project/download/'
path=sorted(glob.glob(path+'*.txt'))

for files in path:

	parse(files)
	printparse()
	#go()
