import string,sys,math,glob,os,sys,re
import numpy as np
np.set_printoptions(threshold=np.nan)
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import protein
from Bio import AlignIO
from scipy import stats
#path='/media/kettina/D874-68C8/alignment/fasta/'
#path='/media/kettina/D874-68C8/alignment/prova/'
#path='/media/D874-68C8/alignment/fasta/'
path='/home/enrichetta/Documents/Disprot/'
path=sorted(glob.glob(path+'*.sth'))

def variable():
	global dis_annotation, ord_annotation, n_dis_region, n_ord_region, av, disar, ordar, hom_used, hom_original, leno, lend, locdis,locord
	dis_annotation={}
	ord_annotation={}
	n_dis_region={}
	n_ord_region={}
	av={}
	disar={}
	ordar={}
	hom_used={}
	hom_original={}
	leno={}
	lend={}
	locdis={}
	locord={}

def entropy_matrix(files):
	
	'''compute the entropy matrix for each sequence starting from multiple sequence alignment'''
	global name,tot,target
	global pos
	global mat, mat_gap
	pos=[]
	drive, path=os.path.splitdrive(files)
     	path, filename=os.path.split(path)
     	name=filename.split('.')[0]
	print name
	align= AlignIO.read(files, 'stockholm')
	tot=len(align)
	hom_original[name]=len(align)
	target=align[0].seq
	alphabet=[]
	alphabet_gap=[]
	pr =set(protein.letters)
	for el in pr:
     		alphabet.append(el)
	alphabet.append('-')
	alphabet_gap.append('-')
		
	for i,el in enumerate(target):
		if el=='-':
			pos.append(i)
	al=[]
	hom_n=0
	for element in align[0:]:
		el=element.seq
		al.append(el)
		hom_n+=1
	hom_used[name]=hom_n

	mat= np.zeros(shape=((len(target)),len(alphabet)))
	mat_gap= np.zeros(shape=((len(target)),len(alphabet_gap)))
	for i in al:
		for x, y in enumerate (i):
			if y in alphabet:
				n=alphabet.index(y)
             			mat[x][n]+=1
				            
         		if y in alphabet_gap:
				n=alphabet_gap.index(y)
				mat_gap[x][n]+=1
			else:
             			pass
			
	for v in range(len(target)):
		mat[v]/=tot
		mat_gap[v]/=tot
	
	for v in range(len(target)):
		g=mat_gap[v][0]
		ind=1-(g)
		mat[v]=(mat[v]*np.log2(mat[v]))/(ind)
		mat[v]=np.absolute(mat[v])
		mat[v]=np.round(mat[v],decimals=5)

def entropy():
	#return a list with the entropies of each sequence
	global entrop	
	sum_h=np.nansum(mat,axis=1)		
	entrop=[]
	for x,y in zip(range(len(sum_h)),sum_h):
		entrop.append(round(y,2))
	offset=0
	for el in pos:
		del entrop[el-offset]
		offset+=1
	
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

def disorder_entropy():
	global loc_dis
	len_d=[]
	loc_dis=[]
	diss=[]
	for element in sorted(dis_regions):	
		n_dis_region[element]=len(dis_regions[element])
		if element==name:
			if  len(dis_regions[element])==0:
				disar[element]=[]
				ln_d=  None
				len_d.append(ln_d)
				lend[element]=len_d
				locdis[element]=loc_dis			
			else:
				for coor in dis_regions[name]:	
					if coor!=' ':
						x=int(coor[0])
						y=int(coor[1])
						ln_d=y-x
						if ln_d>=30:
							len_d.append(ln_d)
							s=0
							n=0
							for i in entrop[x:y+1]: 
								n+=1
								s+=i
							s=s/n
							diss.append(round(s,2))
					else:
						pass
		
				lend[element]=len_d
				disar[element]=diss
				

def order_entropy():	
	ordrr=[]
	len_o=[]
	for element in sorted(ord_regions):
		n_ord_region[element]=len(ord_regions[element])
		if element==name:
			if  len(ord_regions[element])==0:
				ln_o= None
				len_o.append(ln_o)
				leno[name]=len_o
				ordar[element]=ordrr
				
			else:
				for coor in ord_regions[element]:	
					if coor!=' ':
						x=int(coor[0])
						y=int(coor[1])
						ln_o=y-x
						if ln_o >=30:			
							len_o.append(ln_o)
							s=0
							n=0
							for i in entrop[x:y+1]: 
								n+=1
								s+=i
							s=s/n
							ordrr.append(round(s,2)) 
					else:
						pass

				leno[element]=len_o
				ordar[element]=ordrr
				
def normalize_result():
	global minmax
	minmaxdis=[]
	minmaxord=[]
	for k,v in sorted(ordar.iteritems()):
		if len(v)!=0:
			for el in v:
				minmaxord.append(el)
		else:
			pass
	
	for k,v in sorted(disar.iteritems()):
		if len(v)!=0:
			for el in v:
				minmaxdis.append(el)
		else:
			pass
			
	mdis=min(minmaxdis)
	Mdis=max(minmaxdis)
	mord=min(minmaxord)
	Mord=max(minmaxord)
	
	for k,v in sorted(ordar.iteritems()):
		new=[]
		if len(v)!=0:
			for el in v:
				el=round((el-mord)*4.3/float(Mord),2)
				new.append(el)
		ordar[k]=new
	
	for k,v in sorted(disar.iteritems()):
		new=[]
		if len(v)!=0:
			for el in v:
				el=round((el-mdis)*4.3/float(Mdis),2)
				new.append(el)
		disar[k]=new

def average_ord():
	for k,v in sorted(ordar.iteritems()):
		if len(v)==0:
			ord_annotation[k]='*'
		else:
			ord_annotation[k]=round(sum(v)/len(v),2)
	
def average_dis():
	for k,v in sorted(disar.iteritems()):
		if len(v)==0:
			dis_annotation[k]='*'
		else:
			dis_annotation[k]=round(sum(v)/len(v),2)
	
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

def print_result1():
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','generalresult.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:^12} {1:^15} {2:^12} {3:^12} {4:^10}{5:^15}{6:^10}{7:^12}'.format('ID','n_dis_region','n_ord_region','n_hom_orig','n_hom_used','<Hd>','<Ho>','av_entropy')
		for l in sorted(dis_annotation):
			try:
			
				print '{0:^12} {1:^15} {2:^12} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format(l,n_dis_region[l],n_ord_region[l],hom_original[l],hom_used[l],dis_annotation[l],ord_annotation[l],av[l])
			except:
				print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^10}{5:^12}{6:^12}{7:^12}'.format(l ,'*', '*', '*', '*', '*', '*','*')

		p.close()
		sys.stdout=saveout

def print_result2():
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','resultwhitlen.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format('ID','len_dis_reg','dis_H', 'len_ord_reg', 'ord_H','<Hd>','<Ho>','av_entropy')
		for k,v in sorted(leno.iteritems()):
			if k in ordar and lend and disar and dis_annotation and ord_annotation and av:
				for el in map(None,lend[k], disar[k],v, ordar[k]):
					try:
						print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format( k,el[0],el[1], el[2], el[3], dis_annotation[k], ord_annotation[k], av[k])
					except:
						print '{0:^12} {1:^15} {2:^15} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format(name ,'*', '*', '*', '*', '*', '*', '*')
		p.close()
		sys.stdout=saveout

def print_result3():
	
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','disresultlen.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:^12} {1:^15} {2:^15} '.format('ID','len_dis_reg','dis_H')
		for k,v in sorted(lend.iteritems()):
			if k in disar and dis_annotation :
				for el in map(None,v,disar[k]):
					try:
						print '{0:^12} {1:^15} {2:^15} '.format( k,el[0],el[1])
					except:
						print '{0:^12} {1:^15} {2:^15} '.format(k ,'*', '*', )
		p.close()
		sys.stdout=saveout

def print_result4():
	
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','ordresultlen.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:^12} {1:^15} {2:^15} '.format('ID','len_ord_reg' ,'ord_H')
		for k,v in sorted(leno.iteritems()):
			if k in ordar and ord_annotation:
				for el in map(None,v,ordar[k]):
					try:
						print '{0:^12} {1:^15} {2:^15} '.format( k,el[0],el[1])
					except:
						print '{0:^12} {1:^15} {2:^15} '.format(k ,'*', '*', '*')
		p.close()
		sys.stdout=saveout

if __name__ == '__main__':
	regions('/home/enrichetta/Documents/Project/Dataset_oxana/disorder_annotation.fasta')
	#regions('/home/kettina/Scrivania/disorder_annotation.fasta')
	variable()
	for files in path:
		entropy_matrix(files)
		entropy()
		disorder_entropy()
		order_entropy()
	normalize_result()
	average_ord()
	average_dis()
	results()
	print_result1()
	print_result2()
	print_result3()
	print_result4()
