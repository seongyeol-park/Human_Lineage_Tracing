#2019-05-25 add gzip input
import sys, os, gzip
from operator import itemgetter
print(sys.argv[3])
file1=sys.argv[1]
file2=sys.argv[2]
sampleid=sys.argv[3]
caller1=sys.argv[4]
caller2=sys.argv[5]
bin_size=1000000

out_file=open(sampleid+'.'+caller1+'_'+caller2+'_union.vcf','w')
out_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER1\tINFO1\tFORMAT1\tSAMPLE1\tFILTER2\tINFO2\tFORMAT2\tSAMPLE2\tCALLER\n')

chrom_size_dic={}
chr_file=open('/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fa.fai')
chr_line=chr_file.readline().strip()
while chr_line:
	chr_indi=chr_line.split('\t')
	chrom_size_dic[chr_indi[0]]=int(chr_indi[1])
	chr_line=chr_file.readline().strip()
info_list={}
caller_list={}
chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']
print('Only the chromosomes contained in the list below will be reported.')
print(chr_list)

for chrom in chr_list:
	info_list[chrom]={}
	caller_list[chrom]={}
	chrom_size=chrom_size_dic[chrom]
	for i in range(0, chrom_size//bin_size+1):
		info_list[chrom][i]={}
		caller_list[chrom][i]={}

if file1[-3:]=='.gz':
	in_file=gzip.open(file1)
else:
	in_file=open(file1)
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		idx='\t'.join(in_indi[0:2]+in_indi[3:5])
		chr1=in_indi[0]
		pos1=int(in_indi[1])
		info='\t'.join(in_indi[0:10])
		info2='\t'.join(['.','.','.','.'])
		if chr1 not in chr_list:
			in_line=in_file.readline().strip()
			continue
		caller_list[chr1][pos1//bin_size][idx]=caller1
		info_list[chr1][pos1//bin_size][idx]=[info,info2]
	in_line=in_file.readline().strip()
in_file.close()

if file2[-3:] == '.gz':
	in_file=gzip.open(file2)
else:
	in_file=open(file2)
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		idx='\t'.join(in_indi[0:2]+in_indi[3:5])
		chr1=in_indi[0]
		pos1=int(in_indi[1])
		info1='\t'.join(in_indi[0:6])+'\t'+'\t'.join(['.','.','.','.'])
		info2='\t'.join(in_indi[6:10])
		if chr1 not in chr_list:
			in_line=in_file.readline().strip()
			continue
		if idx in caller_list[chr1][pos1//bin_size].keys():
			prev_caller=caller_list[chr1][pos1//bin_size][idx]
			caller_list[chr1][pos1//bin_size][idx]=prev_caller+';'+caller2
			prev_info=info_list[chr1][pos1//bin_size][idx]
			info_list[chr1][pos1//bin_size][idx]=[prev_info[0],info2]
		else:
			caller_list[chr1][pos1//bin_size][idx]=caller2
			info_list[chr1][pos1//bin_size][idx]=[info1,info2]
	in_line=in_file.readline().strip()
in_file.close()

print('Start sorting and writing')

for chrom in chr_list:
	chrom_size=chrom_size_dic[chrom]
	for i in range(0, chrom_size//bin_size +1):
		sorted_list=[]
		for idx in caller_list[chrom][i].keys():
			caller=caller_list[chrom][i][idx]
			info='\t'.join(info_list[chrom][i][idx])
			sorted_list.append([int(idx.split('\t')[1]),info,caller])
		sorted_list.sort(key=itemgetter(0))
		for [pos1, info, caller] in sorted_list:
			out_file.write(info+'\t'+caller+'\n')
