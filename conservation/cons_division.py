import csv
from collections import OrderedDict
def cons(files):
	global consH_Hd,consH_notHd,notconsH_consHd,notconsH_Hd
	f=csv.reader(open(files,'rb'))
	consH_Hd=[]
	consH_notHd=[]
	notconsH_consHd=[]
	notconsH_Hd=[]
	
	for row in f:
		x=row[0].split()[0]
		y=row[0].split()[1]
		z=row[0].split()[2]
		k=row[0].split()[3]
		consH_Hd.append(x)
		consH_notHd.append(y)
		notconsH_consHd.append(z)
		notconsH_Hd.append(k)
	consH_Hd=list(OrderedDict.fromkeys(consH_Hd))
	consH_notHd=list(OrderedDict.fromkeys(consH_notHd))
	notconsH_consHd=list(OrderedDict.fromkeys(notconsH_consHd))
	notconsH_Hd=list(OrderedDict.fromkeys(notconsH_Hd))
	
	del consH_Hd[0], consH_Hd[-1]
	del consH_notHd[0]
	del notconsH_consHd[0], notconsH_consHd[-1]
	del notconsH_Hd[0], notconsH_Hd[-1]
	
def cons2(filex):
	global cons
	global notcons
	f=csv.reader(open(filex,'rb'))
	cons=[]
	notcons=[]
	for row in f:
		x=row[0].split()[0]
		y=row[0].split()[2]
		cons.append(x)
		notcons.append(y)
	cons=list(OrderedDict.fromkeys(cons))
	notcons=list(OrderedDict.fromkeys(notcons))
	del cons[0], cons[-1]
	del notcons[0], notcons[-1]
	
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

def result():
	print '{0:15} {1:12} {2:12} '.format('consH_Hd' ,len(CONconsH_Hd), len(NOTconsH_Hd))
	print '{0:15} {1:12} {2:12} '.format('consH_notHd' , len(CONconsH_notHd), len(NOTconsH_notHd))
	print '{0:15} {1:12} {2:12} '.format('notconsH_consHd', len(CONnotconsH_consHd), len(NOTnotconsH_consHd))
	print '{0:15} {1:12} {2:12} '.format('notconsH_Hd', len(CONnotconsH_Hd), len(NOTnotconsH_Hd))
cons('/home/kettina/Documenti/subset.txt')
cons2('/home/kettina/Documenti/conservation.txt')
analysis()
result()
