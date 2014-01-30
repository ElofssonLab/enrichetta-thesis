from Bio import SeqIO
import re, os, sys
def insert_column(file_name):
	global num_hom	
	num_hom={}
	with open(file_name) as f:
		f=f.readlines()[1:]
		for lines in f:
			lines=lines.split(',')
			num_hom[lines[0]]=lines[7]
	
insert_column('/media/data/dataset_oxana/proteins.txt')

def info_dataset(dis_fasta):
	global disordered
	disordered={}
	global ordered
	ordered={}
	global n_residues
	n_residues={}
	dis_res={}
	ordr_res={}
	unigene={}
	uniprot={}
	n_regions={}
	regions={}
	with open (dis_fasta) as dis_fasta:
		for line in SeqIO.parse(dis_fasta,'fasta'):
	
			info= line.description.split()

			info=info[1].split('|')


			if 'unigene' in info:
				unigene[line.id]= info[info.index('unigene')+1]
			else:
				unigene[line.id]=''	
			if'uniprot' in info:
				
				uniprot[line.id]= info[info.index('uniprot')+1]
			else:
				uniprot[line.id]=''


        		dis=0

        		ordr=0

        		dis+=line.seq.count('0')

        		ordr+=line.seq.count('1')

        		tot=ordr+dis

        		dis_per= round((dis/(float(tot))*100),2)

        		ordr_per=round((ordr/(float(tot))*100),2)

        		dis_res[line.id]=dis

			ordr_res[line.id]=ordr

        		disordered[line.id]=dis_per

        		ordered[line.id]=ordr_per

       			n_residues[line.id]=tot

 			line_str=str(line.seq)
			
			n=1	
			regions[line.id]= [(a.start(),a.end()) for a in list(re.finditer('0+',line_str))]
			
			n_regions[line.id]=len(regions[line.id])
		
			print  '{0:^12} {1:^12} {2:^12} {3:^10} {4:^8} {5:^8} {6:^8} {7:^8} {8:^8} {9:^8} {10:10}'.format(line.id, uniprot[line.id], unigene[line.id], num_hom[line.id],  n_residues[line.id], dis_res[line.id], ordr_res[line.id], disordered[line.id], ordered[line.id], n_regions[line.id], regions[line.id])
			
with open(os.path.join('/home/enrichetta/Project/enrichetta-thesis','dataset.txt'),'w') as f:
	saveout=sys.stdout
        sys.stdout=f
	print '{0:^12} {1:^12} {2:^12} {3:^10} {4:^8} {5:^8} {6:^8} {7:^8} {8:^8} {9:^8} {10:100}'. format('Disprot_ID', 'Uniprot_ID', 'Unigene_ID', 'Num_Hom', 'Tot_Res', 'Dis_Res', 'Ord_Res', '%_Dis', '%_Ord', 'Dis_Reg', 'Coor')
	info_dataset('/media/data/dataset_oxana/disorder_annotation.fasta')
	f.flush()
        f.close()
        sys.stdout=saveout


def general_stat():

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
	print len(key_max_dis)
	for key in key_max_dis:
        	print key, str(max_dis)+'%'
	print
	print 'Max Order'
	print len(key_max_ordr)
	for key in key_max_ordr:
        	print key, str(max_ordr)+'%'

general_stat()



