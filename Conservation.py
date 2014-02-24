import string,sys,math,glob,os,sys,re
import pandas as pd
from pandas import merge
import numpy as np
np.set_printoptions(threshold=np.nan)
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import protein
from Bio import AlignIO
from scipy import stats
path='/media/data/hhblits_files/a2m/'
#path='/home/kettina/Scrivania/Nuova/'
path=sorted(glob.glob(path+'*.a2m'))

dis_annotation={}
ord_annotation={}
n_dis_region={}
n_ord_region={}
av={}


def entropy_matrix(files):
	'''compute the entropy matrix for each sequence starting from multiple sequence alignment'''
	global name
	drive, path=os.path.splitdrive(files)
     	path, filename=os.path.split(path)
     	name=filename.split('.')[0]
	global hom_used
	hom_used={}
	global hom_original
	hom_original={}
	global mat
	align= AlignIO.read(files, 'fasta')
	hom_original[name]=len(align)
	target=align[0].seq
	alphabet=[]
	pr =set(protein.letters)
	for el in pr:
     		alphabet.append(el)

	alphabet.append('-')
	gap=0

	for element in target:
    		gap+=element.count('.')
    		gap+= element.count('-')
	thr=(len(target)-gap)*80/100

	al=[]
	hom_n=0
	for element in align[0:]:
     		gaps=0
     		ins=0
     		ins+=element.seq.count('.')
     		gaps+=element.seq.count('-')
     		inner_thr=len(element.seq)-gaps
     		if inner_thr >=thr:
         		el=element.seq
         		al.append(el)
			hom_n+=1
	hom_used[name]=hom_n
	mat= np.zeros(shape=((len(target)),len(alphabet)))

	for i in al:

		for x, y in enumerate (i):
        


			if y in alphabet:
				n=alphabet.index(y)
             			mat[x][n]+=1
            
         		else:
             			pass


	sumline=np.nansum(mat, axis=1)
	for v in range(len(target)):
		mat[v]/=sumline[v]


	for v in range(len(target)):
		mat[v]=mat[v]*np.log2(mat[v])

		mat[v]=np.absolute(mat[v])
		mat[v]=np.round(mat[v],decimals=5)



	#df=pd.DataFrame(mat, index=[range(len(target))],columns=[alphabet])
	#pd.set_option('display.height',20)
	#pd.set_option('display.max_rows', len(target))
	#pd.set_option('display.max_columns',len(alphabet))
	#pd.set_option('display.precision',3)
	#print df



def entropy():
	#return a list with the entropies of each sequence
	sum_h=np.nansum(mat,axis=1)
	global entrop
	entrop=[]
	for x,y in zip(range(len(sum_h)),sum_h):
	
		entrop.append(round(y,2))
	
	dis_reg=[]
	ord_reg=[]


	for i,j in zip(entrop, range(len(entrop))):
		if math.isnan(float(entrop[j])):
	
			 entrop[j]='#'
		
	
	entrop=[x for x in entrop if x!='#']
	
	

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
                    dis_regions[line.id]= [(a.start(),a.end()) for a in list(re.finditer('0+',line_str))]
                    ord_regions[line.id]= [(a.start(),a.end()) for a in list(re.finditer('1+',line_str))]
        
	#print dis_regions
	#print ord_regions
	    
def disorder_entropy():
	global dis_ar,ln_d,len_d
	len_d=[]
	diss=[]
	for element in sorted(dis_regions):
		
		n_dis_region[element]=len(dis_regions[element])
		if element==name:
			
			if  len(dis_regions[element])==0:
				dis_annotation[name]='*'
				dis_ar=np.array(diss, dtype=np.float64)
				ln_d=  None
				len_d.append(ln_d)
			else:
			
			
				for coor in dis_regions[name]:
					
					if coor!=' ':
						x=int(coor[0])
						y=int(coor[1])
						
						ln_d=y-x
						
							
						len_d.append(ln_d)
							
						s=0
						n=0
						for i in entrop[x:y]: 
							n+=1
							s+=i
						s=s/n
						
						diss.append(round(s,2)) 
						dis_ar=np.array(diss, dtype=np.float64)
	
						l=sum(diss)
						x= len(diss)
						dis_annotation[name]=round((l/x),2)
					else:
						pass
	

def order_entropy():
	global ord_ar,ln_o,len_o				
	ordrr=[]
	len_o=[]
	for element in sorted(ord_regions):
		n_ord_region[element]=len(ord_regions[element])
		if element==name:
			
			if  len(ord_regions[element])==0:
				ord_annotation[element]='*'
				ord_ar=np.array(ordrr, dtype=np.float64)
				ln_o= None
				len_o.append(ln_o)
			else:
				
				for coor in ord_regions[element]:
					
					
					if coor!=' ':
				
						
						x=int(coor[0])
						y=int(coor[1])
					
						ln_o=y-x
							
						len_o.append(ln_o)
							
						s=0
						n=0
						for i in entrop[x:y]: 
							n+=1
							s+=i
						s=s/n
						
						ordrr.append(round(s,2)) 
						ord_ar=np.array(ordrr, dtype=np.float64)
						
						l=sum(ordrr)
						x= len(ordrr)
						ord_annotation[name]=round((l/x),2)
				
					else:
						pass
	
	
	
def t_test():
	
	global t_val, p_val
	
	p_val={}
	t_val={}
	if not ord_ar.size:
		t_val[name]='*'
		p_val[name]='*'
	elif not dis_ar.size:
		t_val[name]='*'
		p_val[name]='*'
	else:
		t,p=stats.ttest_ind(dis_ar,ord_ar,False)
		if np.isnan(t)==True:
			t_val[name]='*'
		else:	
			t_val[name]=round(t,2)
		if np.isnan(p)==True:
			p_val[name]='*'
		else:
			p_val[name]=round(p,2)


def results():
	for k,j in zip(dis_annotation, ord_annotation):
		if dis_annotation[k]!='*' and ord_annotation[j]!='*':
			if dis_annotation[k]==0.0 and ord_annotation[k]==0.0:
				av[k]=0.0
			elif ord_annotation[k]==0.0:
				av[k]=dis_annotation[k]
			else:
				
				av_H=dis_annotation[k]/ord_annotation[k]
				av[k]=round(av_H,2)
			
		elif dis_annotation[k]=='*' and ord_annotation[j]!='*':
			av[k]=ord_annotation[j]
		elif dis_annotation[k]!='*' and ord_annotation[j]=='*':
			av[k]=dis_annotation[k]
					 
	
	#for el in map(None,len_d,dis_ar,len_o,ord_ar):
		#print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format(name ,el[0], el[1], el[2], el[3], dis_annotation[name], ord_annotation[name], av[name])

	#print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}'.format(name,dis_annotation[name],ord_annotation[name], t_val[name], p_val[name])
	print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format(name,n_dis_region[name],n_ord_region[name],hom_original[name],hom_used[name],dis_annotation[name],ord_annotation[name],av[name])

#print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format('ID','len_dis_reg','dis_H', 'len_ord_reg', 'ord_H','<Hd>','<Ho>','av_entropy')
print '{0:^12} {1:^15} {2:^15} {3:^12} {4:^12}{5:^12}{6:^12}{7:^12}'.format('ID','n_dis_region','n_ord_region','n_hom_orig','n_hom_used','<Hd>','<Ho>','av_entropy')
#print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}'.format('ID','<Hd>','<Ho>', 't-value', 'p-value')

regions('/media/data/dataset_oxana/disorder_annotation.fasta')
#regions('/home/kettina/Scrivania/disorder_annotation.fasta')

for files in path:
	
	entropy_matrix(files)
	entropy()
	disorder_entropy()
	order_entropy()
	t_test()
	results()


