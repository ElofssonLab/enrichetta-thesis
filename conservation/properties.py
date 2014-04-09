import string,sys,math,glob,os,sys,re
import numpy as np
np.set_printoptions(threshold=np.nan)
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import protein
from Bio import AlignIO
from scipy import stats
from collections import Counter
#path='/media/kettina/D874-68C8/alignment/fasta/'
#path='/media/kettina/D874-68C8/alignment/prova/'
path='/home/enrichetta/Documents/Disprot/'
#path='/media/data/fasta/'
path=sorted(glob.glob(path+'*.sth'))

def variable():
	global  n_dis_region, n_ord_region, leno, lend, hydrop, pol, posit, neg, specl, gps, disHYDROP, disPOL, disPOSIT, disNEG, disSPECL, disGPS,ordHYDROP, ordPOL, ordPOSIT, ordNEG, ordSPECL, ordGPS
	n_dis_region={}
	n_ord_region={}
	leno={}
	lend={}
	hydrop={}
	pol={}
	posit={}
	neg={}
	specl={}
	gps={}
	disHYDROP={}
	disPOL={}
	disPOSIT={}
	disNEG={}
	disSPECL={}
	disGPS={}
	ordHYDROP={}
	ordPOL={}
	ordPOSIT={}
	ordNEG={}
	ordSPECL={}
	ordGPS={}

def properties_matrix(files):
	
	'''compute the properties matrix for each sequence starting from multiple sequence alignment'''
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
	
	target=align[0].seq
	alphabet=[]
	alphabet_gap=[]
	hydrophobic=['M','A','V','L','I','F','Y','W']
	polar=['S','T','N','Q']
	positive=['K','R','H']
	negative=['D','E']
	special=['P','C','G']
	gap=['-']
	
	for el in hydrophobic:
     		alphabet.append(el)
	for el in polar:
		alphabet.append(el)
	for el in positive:
		alphabet.append(el)
	for el in negative:
		alphabet.append(el)
	for el in special:
		alphabet.append(el)

	alphabet.append('-')
	alphabet_gap.append('-')
	for i,el in enumerate(target):
		if el=='-':
			pos.append(i)

	al=[]
	for element in align[0:]:
     		
         	el=element.seq
         	al.append(el)
	
	mat= np.zeros(shape=((len(target)),len(alphabet)))
	for i in al:
		for x, y in enumerate (i):
			if y in alphabet:
				n=alphabet.index(y)
             			mat[x][n]+=1
			else:
             			pass
	
		
	for v in range(len(target)):
		mat[v]/=tot
		
	n=0
	for v in mat:
		
		hydrop[n]=sum(v[0:7+1])
		pol[n]=sum(v[8:11+1])
		posit[n]=sum(v[12:14+1])
		neg[n]=sum(v[15:16+1])
		specl[n]=sum(v[17:19+1])
		gps[n]=v[20]
		n+=1

def remove():
	#remove the gap position from the different dictionaries containing % of properties
	for el in pos:
		if el in hydrop:
			del hydrop[el]
	
	for el in pos:
		if el in pol:
			del pol[el]
	
	for el in pos:
		if el in posit:
			del posit[el]

	for el in pos:
		if el in neg:
			del neg[el]

	for el in pos:
		if el in specl:
			del specl[el]

	for el in pos:
		if el in gps:
			del gps[el]
	
def updatedict():
	n=0
	for k, v in hydrop.iteritems():
		hydrop[n]=hydrop.pop(k)
		posit[n]=posit.pop(k)
		pol[n]=pol.pop(k)
		neg[n]=neg.pop(k)
		specl[n]=specl.pop(k)
		gps[n]=gps.pop(k)
		n+=1

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

def perdis():
	global loc_dis
	len_d=[]
	loc_dis=[]
	diss=[]
	Hyd=[]
	Pol=[]
	Pos=[]
	Neg=[]
	Spc=[]
	Gp=[]
	for element in sorted(dis_regions):	
		n_dis_region[element]=len(dis_regions[element])
		if element==name:
			if  len(dis_regions[element])==0:
				
				ln_d=  None
				len_d.append(ln_d)
				lend[element]=len_d
				disHYDROP[element]=Hyd
				disPOL[element]=Pol
				disPOSIT[element]=Pos
				disNEG[element]=Neg
				disSPECL[element]=Spc
				disGPS[element]=Gp
					
			else:
				for coor in dis_regions[name]:	
					if coor!=' ':
						x=int(coor[0])
						y=int(coor[1])
						ln_d=y-x
						if ln_d>=30:
							len_d.append(ln_d)
							hyd=[]
							pl=[]
							ps=[]
							ng=[]
							sp=[]
							gp=[]
							
							for k,j in hydrop.iteritems():
								
								if k in range(x,y):
									hyd.append(j)
							try:
								res=round(sum(hyd)/len(hyd)*100,2)
							except:
								res=0
							Hyd.append(res)
							for k,j in pol.iteritems():
								
								if k in range(x,y):
									pl.append(j)
							try:
								res=round(sum(pl)/len(pl)*100,2)
							except:
								res=0
							Pol.append(res)
							for k,j in posit.iteritems():
								
								if k in range(x,y):
									ps.append(j)
							try:
								res=round(sum(ps)/len(ps)*100,2)
							except:
								res=0
							Pos.append(res)
							for k,j in neg.iteritems():
								
								if k in range(x,y):
									ng.append(j)
							try:
								res=round(sum(ng)/len(ng)*100,2)
							except:
								res=0
							Neg.append(res)
							for k,j in specl.iteritems():
								
								if k in range(x,y):
									sp.append(j)
							try:
								res=round(sum(sp)/len(sp)*100,2)
							except:
								res=0
							Spc.append(res)
							for k,j in gps.iteritems():
								
								if k in range(x,y):
									gp.append(j)
							try:
								res=round(sum(gp)/len(gp)*100,2)
							except:
								res=0
							Gp.append(res)
							
				lend[element]=len_d
				disHYDROP[element]=Hyd
				disPOL[element]=Pol
				disPOSIT[element]=Pos
				disNEG[element]=Neg
				disSPECL[element]=Spc
				disGPS[element]=Gp
	
def perord():
	global loc_dis
	len_o=[]
	loc_dis=[]
	diss=[]
	Hyd=[]
	Pol=[]
	Pos=[]
	Neg=[]
	Spc=[]
	Gp=[]
	for element in sorted(ord_regions):	
		n_dis_region[element]=len(ord_regions[element])
		if element==name:
			if  len(ord_regions[element])==0:
				
				ln_o=  None
				len_o.append(ln_o)
				leno[element]=len_o
				
				ordHYDROP[element]=Hyd
				ordPOL[element]=Pol
				ordPOSIT[element]=Pos
				ordNEG[element]=Neg
				ordSPECL[element]=Spc
				ordGPS[element]=Gp
					
			else:
				for coor in ord_regions[name]:	
					if coor!=' ':
						x=int(coor[0])
						y=int(coor[1])
						
						ln_o=y-x
						if ln_o>=30:
							len_o.append(ln_o)
							hyd=[]
							pl=[]
							ps=[]
							ng=[]
							sp=[]
							gp=[]
							for k,j in hydrop.iteritems():
								
								if k in range(x,y):
									
									hyd.append(j)
							
							
							try:
								res=round(sum(hyd)/len(hyd)*100,2)
							except:
								res=0
							Hyd.append(res)
							for k,j in pol.iteritems():
								
								if k in range(x,y):
									pl.append(j)
							try:
								res=round(sum(pl)/len(pl)*100,2)
							except:
								res=0
							Pol.append(res)
							for k,j in posit.iteritems():
								
								if k in range(x,y):
									ps.append(j)
							try:
								res=round(sum(ps)/len(ps)*100,2)
							except:
								res=0
							Pos.append(res)
							for k,j in neg.iteritems():
								
								if k in range(x,y):
									ng.append(j)
							try:
								res=round(sum(ng)/len(ng)*100,2)
							except:
								res=0
							Neg.append(res)
							for k,j in specl.iteritems():
								
								if k in range(x,y):
									sp.append(j)
							try:
								res=round(sum(sp)/len(sp)*100,2)
							except:
								res=0
							Spc.append(res)
							for k,j in gps.iteritems():
								
								if k in range(x,y):
									gp.append(j)
							try:
								res=round(sum(gp)/len(gp)*100,2)
							except:
								res=0
							Gp.append(res)
							
				leno[element]=len_o
				ordHYDROP[element]=Hyd
				ordPOL[element]=Pol
				ordPOSIT[element]=Pos
				ordNEG[element]=Neg
				ordSPECL[element]=Spc
				ordGPS[element]=Gp

def print_result1():
	
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','disresultproperty.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:^12} {1:^12} {2:^12} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format('ID','len_dis','%Hydrop', '%polar','%pos_ch', '%neg_ch','%Special','%gaps')
		for k,v in sorted(lend.iteritems()):
			
			if k in  disHYDROP and disPOL and disPOSIT and disNEG and disSPECL and disGPS and ordHYDROP and ordPOL and ordPOSIT and ordNEG and ordSPECL and ordGPS:
				for el in map(None,v, disHYDROP[k],disPOL[k],disPOSIT[k],disNEG[k],disSPECL[k],disGPS[k]):
					try:
						print '{0:^12} {1:^12} {2:^12} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format( k,el[0],el[1], el[2], el[3],el[4], el[5],el[6])
					except:
						print '{0:^12} {1:^12} {2:^12} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format(name ,'*', '*', '*', '*', '*', '*', '*')
		p.close()
		sys.stdout=saveout
def print_result2():
	
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','ordresultproperty.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:^12} {1:^12} {2:^12} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format('ID','len_ord','%Hydrop', '%polar','%pos_ch', '%neg_ch','%Special','%gaps')
		for k,v in sorted(leno.iteritems()):
			
			if k in ordHYDROP and ordPOL and ordPOSIT and ordNEG and ordSPECL and ordGPS:
				for el in map(None,v, ordHYDROP[k],ordPOL[k],ordPOSIT[k], ordNEG[k], ordSPECL[k],ordGPS[k]):
					try:
						print '{0:^12} {1:^12} {2:^12} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format( k,el[0],el[1], el[2], el[3],el[4], el[5],el[6])
					except:
						print '{0:^12} {1:^12} {2:^12} {3:^12}{4:^12}{5:^12}{6:^12}{7:^12}'.format(name ,'*', '*', '*', '*', '*', '*', '*')
		p.close()
		sys.stdout=saveout

if __name__ == '__main__':
	regions('/home/enrichetta/Documents/Project/Dataset_oxana/disorder_annotation.fasta')
#	regions('/home/kettina/Scrivania/disorder_annotation.fasta')
	variable()
	for files in path:
		properties_matrix(files)	
		remove()
		updatedict()
		perdis()
		perord()
	print_result1()
	print_result2()		
