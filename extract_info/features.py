from Bio import SwissProt
import glob
from Bio import SeqIO
import re
import csv
import urllib2
import os,sys
import rdflib

def ids(files):
	f=csv.reader(open(files,'rb'))
	global dic_id
	dic_id={}
	for row in f:
		x=row[0].split()
		dic_id[x[0]]=x[1]
	del dic_id['Disprot_ID']

def extract_col(files):
	f=csv.reader(open(files,'rb'))
	for row in f:
		x=row[0].split()	
	col=len(x)
	
	return col
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
					with open(os.path.join('/home/kettina/Scrivania/Project/Download',file), 'w') as f:
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
								with open(os.path.join('/home/kettina/Scrivania/Project/Download',d+'.txt'), 'w') as f:
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
			
	with open(os.path.join('/home/kettina/Scrivania/Project/Results/Comparison/'+file_name,file_name+'_nofile.txt'),'w') as d:
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
		with open(os.path.join('/home/kettina/Scrivania/Project/Results/Comparison/'+file_name,file_name+'_nofile.txt'),'a') as d:
			saveout=sys.stdout
        		sys.stdout=d
			print em_f
			d.flush()
        		d.close()
        		sys.stdout=saveout
	else:
		pass
	for el in em_f:
		print '/home/kettina/Scrivania/Project/Download/'+el+'.txt'
		os.remove('/home/kettina/Scrivania/Project/Download/'+el+'.txt')	

def regions(dis_fasta):
	#extract disordered/ordered regions from fasta file
	global dis_regions
	global ord_regions
	dis_regions={}
        ord_regions={}
        with open (dis_fasta) as dis_fasta:
                for line in SeqIO.parse(dis_fasta,'fasta'):
                    info= line.description.split()
                    info=info[1].split('|')
		    
                    line_str=str(line.seq)
                    dis_regions[info[3]]= [(a.start(),a.end()) for a in list(re.finditer('0+',line_str))]
                    ord_regions[info[3]]= [(a.start(),a.end()) for a in list(re.finditer('1+',line_str))]
	
	#for k, v in dis_regions.iteritems():
	#	for coor in v:
	#		print k, dic_id[k],coor[0], coor[1]
def features(files):
	ft=['ZN_FING', 'REGION','METAL','SITE','SIGNAL','REPEAT', 'NP_REGION', 'BINDING','MOTIF','MOD_RES', 'LIPID','DOMAIN','DNA_BIND','DISULFID','CROSSLNK', 'CARBOHYD','CA_BIND', 'ACT_SITE']
	for record in SwissProt.parse(open(files)):
		for l in record.features:
			
			if l[0] in ft:
				print l[0]+','+str(l[1])+'-'+str(l[2])+','+l[3]


def dic_coordis(files):
	global dicdis
	dicdis={}
	coor={}
	f=csv.reader(open(files,'rb'))
	for k,v in dic_id.iteritems():
		dicdis[k]=coor
	
	
	k1=[]
	k2=[]
	v=[]
	
	for row in f:
		
		k1.append(row[0].split()[0])

		k2.append(row[0].split()[1])
		v.append(row[0].split()[3])
	for n in enumerate(k):
		print n
	'''		
		if k in dicdis:
			dicdis[k][kb]='x'
	print dicdis

		if dicdis.has_key(k1):
			dicdis[k1][k2]=v
	print dicdis

		dicdis[row[0].split()[0]]={}
		if row[0].split()[0] in dicdis:
			dicdis[row[0].split()[0]][row[0].split()[1]]='x'
	print dicdis
	
		for k in dicdis.keys():
			if row[0].split()[0]== k:
				print row[0].split()[0], k, row[0].split()[1]
				dicdis[k][row[0].split()[1]]='x'#row[0].split()[3]
	#print dicdis
	
		#print row[0].split()[3]
		#dicdis[row[0].split()[0]]={}
			dicdis[row[0].split()[0]][row[0].split()[1]]=row[0].split()[3]
	del dicdis['ID']
	print dicdis
	'''
def parse_features_dis(files):
	f=csv.reader(open(files,'rb'))
	drive, path=os.path.splitdrive(files)
    	path, filename=os.path.split(path)
	c=filename.split('.')[0]
	
	for row in f:
		x=str(row[1].split('-')[0])
		y=row[1].split('-')[1]
		if '?' in x:
			x=x.replace('?','')
		
		if '?' in y:
			y=y.replace('?','')
			
		for k, v in dis_regions.iteritems():
		
			if c==k:
				for coor in v:
					coor1=int(coor[0])
					coor2=int(coor[1])
					lncor=coor2-coor1
					if lncor >=30:
				#print k,dis_regions[k], coor
				#print k, dic_id[k],coor[0]+'-'+ coor[1]
						if int(x) >= coor[0] and int(y) <= coor[1]:
							xx=dic_id.keys()[dic_id.values().index(k)].split()
							
							for l in xx:
								print l+','+k+','+str(x)+'-'+str(y)+','+str(coor[0])+'-'+str(coor[1])+','+ row[0]+','+ row[2]
def add_entropy_dis(filex, files):
	f=csv.reader(open(files,'rb'))
	d=csv.reader(open(filex,'rb'))
	k1=[]
	k2=[]
	k3=[]
	for row in d:
		row=row[0].split()
		k1.append(row[0])
		k2.append(row[1])
		k3.append(row[3])
	for row in f:
		for k,v,j in zip(k1,k2,k3):
			if row[0]==k:
				if row[3]==v:
					print row[0]+','+row[1]+','+row[2]+','+row[3]+','+row[4]+','+row[5]+','+j

def parse_features_ord(files):
	f=csv.reader(open(files,'rb'))
	drive, path=os.path.splitdrive(files)
    	path, filename=os.path.split(path)
	c=filename.split('.')[0]
	
	for row in f:
		x=str(row[1].split('-')[0])
		y=row[1].split('-')[1]
		if '?' in x:
			x=x.replace('?','')
		
		if '?' in y:
			y=y.replace('?','')
			
		for k, v in ord_regions.iteritems():
		
			if c==k:
				for coor in v:
					coor1=int(coor[0])
					coor2=int(coor[1])
					lncor=coor2-coor1
					if lncor >=30:
				#print k,dis_regions[k], coor
				#print k, dic_id[k],coor[0]+'-'+ coor[1]
						if int(x) >= coor[0] and int(y) <= coor[1]:
							xx=dic_id.keys()[dic_id.values().index(k)].split()
							
							for l in xx:
								print l+','+k+','+str(x)+'-'+str(y)+','+str(coor[0])+'-'+str(coor[1])+','+ row[0]+','+ row[2]

def add_entropy_ord(filex, files):
	f=csv.reader(open(files,'rb'))
	d=csv.reader(open(filex,'rb'))
	k1=[]
	k2=[]
	k3=[]
	for row in d:
		row=row[0].split()
		k1.append(row[0])
		k2.append(row[1])
		k3.append(row[3])
	for row in f:
		for k,v,j in zip(k1,k2,k3):
			if row[0]==k:
				if row[3]==v:
					print row[0]+','+row[1]+','+row[2]+','+row[3]+','+row[4]+','+row[5]+','+j
		
if __name__=='__main__':
	ids('/home/kettina/Documenti/enrichetta-thesis/dataset.txt')
	col=extract_col('/home/kettina/Scrivania/Project/Results/subsetHd.txt')
	
	regions('/home/kettina/Scrivania/disorder_annotation.fasta')
	
	for n in range(col):
		info('/home/kettina/Scrivania/Project/Results/subsetHd.txt',n) #select the column number
		
		#download()
		
	path='/home/kettina/Scrivania/Project/Download/'
	path=sorted(glob.glob(path+'*.txt'))
	#print path
	
	for files in path:
		
		empty(files)
		drive, path=os.path.splitdrive(files)
    		path, filename=os.path.split(path)
		c=filename
		#print c
		with open(os.path.join('/home/kettina/Scrivania/Project/Download2/',c),'w') as f:
			saveout=sys.stdout
        		sys.stdout=f
			features(files)
			f.flush()
        		f.close()
        		sys.stdout=saveout
	path2='/home/kettina/Scrivania/Project/Download2/'
	path2=sorted(glob.glob(path2+'*.txt'))	
	#dic_coordis('/home/kettina/Scrivania/prova/disresultlen.txt')
	with open(os.path.join('/home/kettina/Scrivania/prova','disorderfeatures.txt'),'w') as f:
			saveout=sys.stdout
        		sys.stdout=f
			for files in path2:
		
				parse_features_dis(files)
	
			
			f.flush()
        		f.close()
        		sys.stdout=saveout
	with open(os.path.join('/home/kettina/Scrivania/prova','disorderfeatures+entropy.txt'),'w') as f:
		saveout=sys.stdout
        	sys.stdout=f
		add_entropy_dis('/home/kettina/Scrivania/prova/disresultlen.txt','/home/kettina/Scrivania/prova/disorderfeatures.txt')
		
		f.flush()
        	f.close()
       		sys.stdout=saveout
	with open(os.path.join('/home/kettina/Scrivania/prova','orderfeatures.txt'),'w') as f:
			saveout=sys.stdout
        		sys.stdout=f
			for files in path2:
		
				parse_features_ord(files)
	
			
			f.flush()
        		f.close()
        		sys.stdout=saveout
	with open(os.path.join('/home/kettina/Scrivania/prova','orderfeatures+entropy.txt'),'w') as f:
		saveout=sys.stdout
        	sys.stdout=f
		add_entropy_ord('/home/kettina/Scrivania/prova/ordresultlen.txt','/home/kettina/Scrivania/prova/orderfeatures.txt')
		
		f.flush()
        	f.close()
       		sys.stdout=saveout
