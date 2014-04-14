import csv
import urllib2
import os,sys
import glob
from  goatools import obo_parser 
import rdflib
from Bio import SwissProt

def ids(files):
	f=csv.reader(open(files,'rb'))
	global dic_id
	dic_id={}
	for row in f:
		x=row[0].split()
		dic_id[x[0]]=x[1]
	del dic_id['Disprot_ID']
	
def info(files,n): # n is the column number that you want extract
	global ids,file_name,temp_id
	ids=[]	
	f=csv.reader(open(files,'rb'))
	for row in f:
		x=row[0].split()
		
		if x[n]!='*':
			ids.append(x[n])
	file_name=ids[0]
	print file_name
	
	del ids[0]
	temp_id={k:dic_id[k] for k in ids if k in dic_id}
	
def download():
	no_file=[]
	new=[]
	global el
	for el in ids:
		if el in temp_id.keys():
			if dic_id[el] !='N/A':	
				#print '{0:12}'.format(temp_id[el]), el
				x=[]
				file=dic_id[el]+'.txt'
				fil=dic_id[el]+'.rdf'
				url='http://www.uniprot.org/uniprot/'+temp_id[el]+'.txt'
			
				try:
					response=urllib2.urlopen(url)
					#print url
					with open(os.path.join('/home/enrichetta/Documents/Project/Download',file), 'w') as f:
						f.write(response.read())
				except:
					no_file.append(temp_id[el])
					print no_file
					try:
						g=rdflib.Graph()
						g.load('http://www.uniprot.org/uniprot/'+temp_id[el]+'.rdf')
						l=[]
						for s,p,o in g:
							l.append(o)
						for z in l:	
							
							
							if '/uniprot/' in z:
								d=z.split('/')[4]
								x.append(d)
								
								url='http://www.uniprot.org/uniprot/'+d+'.txt'
								response=urllib2.urlopen(url)
								with open(os.path.join('/home/enrichetta/Documents/Project/Download',d+'.txt'), 'w') as f:
									f.write(response.read())
						x=[str(l) for l in x]
						print x
						for k,v in temp_id.iteritems():
							if v in no_file:
								temp_id[k]=x
							else:
								pass
						for n in x:
							new.append(n)
					except:
						pass
					
					pass
				#with open(os.path.join('/home/enrichetta/Project/download',file), 'w') as f:
				#	f.write(response.read())
			else:
				no_file.append(temp_id[el])
		else:
			no_file.append(temp_id[el])
			
	with open(os.path.join('/home/enrichetta/Documents/Project/Results/Comparison/'+file_name,file_name+'nofile.txt'),'w') as d:
		saveout=sys.stdout
        	sys.stdout=d
		
		for l in no_file:
			
			print l
		d.flush()
        	d.close()
        	sys.stdout=saveout	
		
def empty(files):
	global em_f
	em_f=[]
	drive, path=os.path.splitdrive(files)
     	path, filename=os.path.split(path)
     	name=filename.split('.')[0]
	if os.stat(files).st_size == 0:
    		em_f.append(name)
		with open(os.path.join('/home/enrichetta/Documents/Project/Results/Comparison/'+file_name,file_name+'nofile.txt'),'a') as d:
			saveout=sys.stdout
        		sys.stdout=d
			print em_f
			d.flush()
        		d.close()
        		sys.stdout=saveout
	else:
		pass
	for el in em_f:
		print '/home/enrichetta/Documents/Project/Download/'+el+'.txt'
		os.remove('/home/enrichetta/Documents/Project/Download/'+el+'.txt')
def parse(files):
	global name,spec
	global GOs,Pfam,Ensembl
	global species,spec,tax
	global idss
	spec={}
	drive, path=os.path.splitdrive(files)
     	path, filename=os.path.split(path)
     	name=filename.split('.')[0]
	gOs=[]
	GOs={}
	Pfam={}
	Ensembl={}
	pfam=[]
	ensembl=[]
	tax={}
	species={}
	idss={}
	x=''
	try:
		
		for record in SwissProt.parse(open(files)):
			species[name]=record.organism[0:-1]
			idss[name]=record.entry_name
			spec[name]=record.organism_classification[-1]
			tax[name]=record.organism_classification[0]
			
			x=len(record.cross_references)
			c=0
			d=0
			
    			for l in record.cross_references:
				
				if l[0]=='Pfam':
					
					c+=1
					if x-c!=x:
						pfam.append(l[1])
				if x-c==x:	
					pfam.append('*')
				if l[0]=='GO':
					gOs.append(l[1])
				if l[0]=='Ensembl':
					
					d+=1
					if x-d!=x:
						ensembl.append(l[3])
				
				
				
				
				if l[0]=='EnsemblBacteria':
				
					d+=1
					if x-d!=x:
						ensembl.append(l[3])
				if l[0]=='EnsemblFungi':
					d+=1
					if x-d!=x:
						ensembl.append(l[3])
				
			
				if l[0]=='EnsemblMetazoa':
					d+=1
					if x-d!=x:
						ensembl.append(l[3])
				
				if l[0]=='EnsemblPlants':
					d+=1
					if x-d!=x:
						ensembl.append(l[3])
				
				if l[0]=='EnsemblProtists':
					d+=1
					if x-d!=x:
						ensembl.append(l[3])
				if x-d==x:
					ensembl.append('*')
				
				
				
		x=len(sorted(set(pfam)))
		
		pfam=sorted(set(pfam))
		ensembl=sorted(set(ensembl))
		
		if x!=1:
			del pfam[0]
		for l in ensembl:
			if len(ensembl)!= 1 and l=='*':
				del ensembl[0]
		Pfam[name]=pfam
		Ensembl[name]=ensembl
		
		GOs[name]=gOs
		'''
		for l in EnsemblBacteria:
			if len(EnsemblBacteria)!= 1 and l=='*':
				del EnsemblBacteria[0]
		
		for l in EnsemblFungi:
			if len(EnsemblFungi)!= 1 and l=='*':
				del EnsemblFungi[0]
		for l in EnsemblMetazoa:
			if len(EnsemblMetazoa)!= 1 and l=='*':
				del EnsemblMetazoa[0]
		for l in EnsemblPlants:
			if len(EnsemblPlants)!= 1 and l=='*':
				del (EnsemblPlants)[0]
		for l in EnsemblProtists:
			if len (EnsemblProtists)!=1 and l=='*':
				del EnsemblProtists[0]
		'''
	except:
		pass

	
def printpfam():
	tmp_dic={}
	for key,value in temp_id.items():
		if value==name:
			
			if value in Pfam:
				for el in Pfam[value]:
					print'{0:12} {1:12} {2:12} '.format( key,value, el)
		if type(value) is list:
			for j in value:
				tmp_dic[j]=key
	
	for key,value in sorted(tmp_dic.items()):
			
			if key in Pfam:
				
				for l in Pfam[key]:

				
					print'{0:12} {1:12} {2:12} '.format( value,key,l)
		
def printpfam2():
	tmp_dic={}
	for key,value in temp_id.items():
		if type(value) is list:
			for j in value:
				tmp_dic[j]=key
	
	for key,value in sorted(tmp_dic.items()):
			
			if key in Pfam:
				
				for l in Pfam[key]:

				
					print'{0:12} {1:12} {2:12} '.format( value,key,l)
			
def printgo():	
	global sp,p
	tmp_dic={}			
	#sp={'cellular_component':'C', 'biological_process':'P', 'molecular_function':'F'}
	#p=obo_parser.GODag('/home/enrichetta/Documents/goatools/data/gene_ontology.1_2.obo')
	
	for key,value in temp_id.items():
		if value==name:
			if value in GOs:
				for l in GOs[value]:
		
					print '{0:12} {1:12} {2:15} '.format(key,value,l)#, sp[p[l].namespace], p[l].name)
		if type(value) is list:
			for j in value:
				tmp_dic[j]=key
	for key,value in sorted(tmp_dic.items()):
		
			if key in GOs:
				for l in GOs[key]:
		
					print '{0:12} {1:12} {2:15} '.format(value,key,l)#, sp[p[l].namespace], p[l].name)

def printspecies():
	tmp_dic={}
	for key,value in temp_id.items():
		if value==name:
			if value in species:
				for val in species:
					#print'{0:12} {1:12} {2:15} {3:15}{4:25} {5:20}'.format( key,value,idss[value],tax[value],spec[value],species[value])
					print key+','+value+','+idss[value]+','+tax[value]+','+spec[value]+','+species[value]
		if type(value) is list:
			for j in value:
				tmp_dic[j]=key
	for key,value in sorted(tmp_dic.items()):
		
			if key in species:
				for l in species:#[key]:
					#print'{0:12} {1:12} {2:15} {3:15}{4:25} {5:20}'.format(value,key,idss[key],tax[key],spec[key],species[key])
					print value+','+key+','+idss[key]+','+tax[key]+','+spec[key]+',',species[key]
def printensembl():
	tmp_dic={}
	for key,value in temp_id.items():
		if value==name:
			if value in Ensembl:
				for l in Ensembl[value]:
				
					print'{0:12} {1:12} {2:15} '.format( key,value, l)
		if type(value) is list:
			for j in value:
				tmp_dic[j]=key
	for key,value in sorted(tmp_dic.items()):
		
			if key in Ensembl:
				for l in Ensembl[key]:
					print'{0:12} {1:12} {2:15} '.format(value,key, l)
		
def delete():
	fil=glob.glob('/home/enrichetta/Documents/Project/Download/*')
	fil2=glob.glob('/home/enrichetta/Documents/Project/Download2/*')
	for f in fil:
		os.remove(f)
	for f in fil2:
		os.remove(f)
		

if __name__=='__main__':
	ids('/home/enrichetta/Documents/Project/Results/dataset.txt')
	
	for n in range(4):
		info('/home/enrichetta/Documents/Project/Results/subset.txt',n) #select the column number
		
		download()
		
		path='/home/enrichetta/Documents/Project/Download/'
		path=sorted(glob.glob(path+'*.txt'))
		path2='/home/enrichetta/Documents/Project/Download2/'
		path2=sorted(glob.glob(path2+'*.txt'))
		
		for files in path:
			empty(files)
			
			
		
		with open(os.path.join('/home/enrichetta/Documents/Project/Results/Comparison/'+file_name,file_name+'GO.txt'),'w') as f:
			saveout=sys.stdout
        		sys.stdout=f
			print '{0:12} {1:12} {2:15} '.format('ID','Uniprot', 'GO_term')#, 'sp','description')
			for files in path:

				parse(files)
				printgo()
			
			f.flush()
        		f.close()
        		sys.stdout=saveout
		
		with open(os.path.join('/home/enrichetta/Documents/Project/Results/Comparison/'+file_name,file_name+'Pfam.txt'),'w') as f:
			saveout=sys.stdout
        		sys.stdout=f
			print'{0:12} {1:12} {2:12} '.format('ID','Uniprot', 'Pfam')
			
			for files in path:

				parse(files)

				printpfam()
			
			f.flush()
        		f.close()
        		sys.stdout=saveout
		
		with open(os.path.join('/home/enrichetta/Documents/Project/Results/Comparison/'+file_name,file_name+'species.txt'),'w') as f:
			saveout=sys.stdout
        		sys.stdout=f
			#print'{0:12} {1:12} {2:15} {3:15}{4:25} {5:20}'.format('ID','Uniprot','idss','Domain','Genus','species')
			print 'ID'+','+'Uniprot'+','+'idss'+','+'Domain'+','+'Genus'+','+'species'
			for files in path:

				parse(files)
			
				printspecies()
		
			f.flush()
        		f.close()
        		sys.stdout=saveout
		with open(os.path.join('/home/enrichetta/Documents/Project/Results/Comparison/'+file_name,file_name+'ensembl.txt'),'w') as f:
			saveout=sys.stdout
        		sys.stdout=f
			print'{0:12} {1:12} {2:15} '.format('ID','Uniprot','Ensembl')
			for files in path:

				parse(files)
			
				printensembl()
		
			f.flush()
        		f.close()
        		sys.stdout=saveout
		delete()
		
	
