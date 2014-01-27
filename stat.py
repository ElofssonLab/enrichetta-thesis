from Bio import SeqIO
dis_fasta=open('/media/data/dataset_oxana/disorder_annotation.fasta','rU')
disordered={}
ordered={}
n_residues={}
for line in SeqIO.parse(dis_fasta,'fasta'):
	dis=0
	ordr=0
	print line.id
		
	dis+=line.seq.count('1')

	print 'disordered residues:', dis
	ordr+=line.seq.count('0')
	print 'ordered residues:',ordr
	tot=ordr+dis
	print 'total residues:',tot
	dis_per= round((dis/(float(tot))*100),2)
	ordr_per=round((ordr/(float(tot))*100),2)
	print 'disordered', dis_per, '%'
	print 'ordered', ordr_per, '%'
	disordered[line.id]=dis_per
	ordered[line.id]=ordr_per
	n_residues[line.id]=tot

tot_proteins=len(disordered)
print
print 'Total number of proteins:', tot_proteins
max_length=max(n_residues.values())
key_max_len=[x for x,y in n_residues.items() if y==max_length]

min_length=min(n_residues.values())
key_min_len=[x for x,y in n_residues.items() if y==min_length]

av_len=sum(n_residues.values())/tot_proteins

max_dis=max(disordered.values())
key_max_dis=[x for x,y in disordered.items() if y==max_dis]

max_ordr=max(ordered.values())
key_max_ordr=[x for x,y in ordered.items() if y==max_ordr]

print
print 'Max length'
print

for key in key_max_len:
	print key, max_length
print
print 'Min length'
print
for key in key_min_len:
	print key, min_length
print
print 'Average length'
print
print av_len
print 
print 'Max Disorder'
for key in key_max_dis:
	print key, str(max_dis)+'%'
print
print 'Max Order'
for key in key_max_ordr:
	print key, str(max_ordr)+'%'

