import csv
import os,sys
import matplotlib.pyplot as plt
import itertools
def info(files):
	global tot, cons, notcons,res, dic_entrop,Hd
	global consH_Hd,consH_notHd,notconsH_consHd,notconsH_Hd,cons_Hd,notcons_Hd
	dic_entrop={}
	dic_Hd={}
	dic_H={}
	dic_region={}
	consH_Hd=[]
	notconsH_consHd=[]
	consH_notHd=[]
	notconsH_Hd=[]
	cons=[]
	notcons=[]
	cons_Hd=[]
	notcons_Hd=[]
	f=csv.reader(open(files,'rb'))
	for row in f:
		x=row[0].split()
		
		dic_Hd[x[0]]=x[5]
		dic_H[x[0]]=x[6]
		dic_entrop[x[0]]=x[7]
	
	
	del dic_Hd['ID']
	del dic_H['ID']
	del dic_entrop['ID']
	for x in dic_H:
		if dic_H[x]!='*':
			dic_H[x]=float(dic_H[x])
	for y in dic_Hd:
		
		if dic_Hd[y]!='*':
			
			dic_Hd[y]=float(dic_Hd[y])
	for z in dic_entrop:
		if dic_entrop[z]!='*':
			dic_entrop[z]=float(dic_entrop[z])
	
	
	for k,v in dic_entrop.iteritems():
		if v <=2:
			cons.append(k)
		else:
			notcons.append(k)
	for x, y in zip(dic_H, dic_Hd):
		if dic_H[x]<=2 and dic_Hd[y]<=2:
			consH_Hd.append(y)	
		if dic_H[x]<=2 and dic_Hd[y]>2:
			consH_notHd.append(y)
		if dic_H[x]>2 and dic_Hd[y]<=2:
			notconsH_consHd.append(y)
		if dic_H[x]>2 and dic_Hd[y]>2:
			notconsH_Hd.append(y)
	for k,v in dic_Hd.iteritems():
		if v <=2:
			cons_Hd.append(k)
		else:
			notcons_Hd.append(k)

	tot=itertools.izip_longest(consH_Hd,consH_notHd,notconsH_consHd,notconsH_Hd,fillvalue='*')
	res=itertools.izip_longest(cons, notcons,fillvalue='*')
	Hd=itertools.izip_longest(cons_Hd, notcons_Hd,fillvalue='*')
def print_result1():
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','subset.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:20} {1:20} {2:20} {3:20}'.format('consH_Hd','consH_notHd','notconsH_consHd','notconsH_Hd')
		for el in tot:
		
			print '{0:20} {1:20} {2:20} {3:20}'.format(el[0],el[1],el[2],el[3])
		p.close()
		sys.stdout=saveout

def print_result2():
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','subsetHd.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:20} {1:20}'.format('cons_Hd','notcons_Hd')
		for el in Hd:
		
			print '{0:20} {1:20} '.format(el[0],el[1])
		p.close()
		sys.stdout=saveout
	
def print_result3():
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','conservation.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print'{0:^30} {1:^35}'.format('Conserved','Not Conserved')
		print '{0:^15} {1:^15} {2:^15} {3:^15}'.format('ID','<H>','ID','<H>')
		for el in res:
			if el[0]=='*' and el[1]=='*':
				print '{0:^15} {1:^15} {2:^15} {3:^15}'.format(el[0],'*',el[1],'*')
			elif el[0]=='*' and el[1]!='*':
				print '{0:^15} {1:^15} {2:^15} {3:^15}'.format(el[0],'*',el[1],dic_entrop[el[1]])
			elif el[0]!='*' and el[1]=='*':
				print '{0:^15} {1:^15} {2:^15} {3:^15}'.format(el[0],dic_entrop[el[0]],el[1],'*')
			else:
				print '{0:^15} {1:^15} {2:^15} {3:^15}'.format(el[0],dic_entrop[el[0]],el[1],dic_entrop[el[1]])
			
		p.close()
		sys.stdout=saveout

	
def analysis():
	global CONconsH_Hd, NOTconsH_Hd, CONconsH_notHd, NOTconsH_notHd, CONnotconsH_consHd, NOTnotconsH_consHd, CONnotconsH_Hd, NOTnotconsH_Hd
	CONconsH_Hd=[]
	NOTconsH_Hd=[]
	CONconsH_notHd=[]
	NOTconsH_notHd=[]
	CONnotconsH_consHd=[]
	NOTnotconsH_consHd=[]
	CONnotconsH_Hd=[]
	NOTnotconsH_Hd=[]
	for el in consH_Hd:
		if el in cons:
			CONconsH_Hd.append(el)
		if el in notcons:
			NOTconsH_Hd.append(el)
	for el in consH_notHd:
		if el in cons:
			CONconsH_notHd.append(el)
		if el in notcons:
			NOTconsH_notHd.append(el)
	for el in notconsH_consHd:
		if el in cons:
			CONnotconsH_consHd.append(el)
		if el in notcons:
			NOTnotconsH_consHd.append(el)
	for el in notconsH_Hd:
		if el in cons:
			CONnotconsH_Hd.append(el)
		if el in notcons:
			NOTnotconsH_Hd.append(el)
	print CONconsH_notHd
def analysis2():
	global CONcons_Hd, NOTcons_Hd,NOTnotcons_Hd,CONnotcons_Hd
	CONcons_Hd=[]
	NOTcons_Hd=[]
	CONnotcons_Hd=[]
	NOTnotcons_Hd=[]
	for el in cons_Hd:
		if el in cons:
			CONcons_Hd.append(el)
		if el in notcons:
			NOTcons_Hd.append(el)
	for el in notcons_Hd:
		if el in cons:
			CONnotcons_Hd.append(el)
		if el in notcons:
			NOTnotcons_Hd.append(el)
	
def result():
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','Consnotcons.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:15} {1:12} {2:12} '.format('consH_Hd' ,len(CONconsH_Hd), len(NOTconsH_Hd))
		print '{0:15} {1:12} {2:12} '.format('consH_notHd' , len(CONconsH_notHd), len(NOTconsH_notHd))
		print '{0:15} {1:12} {2:12} '.format('notconsH_consHd', len(CONnotconsH_consHd), len(NOTnotconsH_consHd))
		print '{0:15} {1:12} {2:12} '.format('notconsH_Hd', len(CONnotconsH_Hd), len(NOTnotconsH_Hd))	
		p.close()
		sys.stdout=saveout

def result2():
	with open(os.path.join('/home/enrichetta/Documents/Project/Results','ConsnotconsHd.txt'),'w') as p:
		saveout=sys.stdout
		sys.stdout=p
		print '{0:15} {1:12} {2:12} '.format('cons_Hd' ,len(CONcons_Hd), len(NOTcons_Hd))
		
		print '{0:15} {1:12} {2:12} '.format('notcons_Hd', len(CONnotcons_Hd), len(NOTnotcons_Hd))
		
		p.close()
		sys.stdout=saveout

if __name__ == '__main__':
	info('/home/enrichetta/Documents/Project/Results/generalresult.txt')
	print_result1()	
	print_result2()
	analysis()
	analysis2()
	result()
	result2()
