#2018-06-02 Specification option for sorting
#2019-12-09 write only univerally shared variants
#2020-03-08 if chr1 not in chr_list -> continue
import sys, os
from operator import itemgetter

f_list=open(sys.argv[1])   # path of input files
id_list=open(sys.argv[2])  # id of input files
out_file=open(sys.argv[3],'w') # name of output file
InCaseNum=int(sys.argv[4]) # number of input cases

chrom_size_dic={}
chr_file=open('/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fa.fai')
chr_line=chr_file.readline().strip()
while chr_line:
	chr_indi=chr_line.split('\t')
	chrom_size_dic[chr_indi[0]]=int(chr_indi[1])
	chr_line=chr_file.readline().strip()
row_list={}
case_list={}
bin_size=500000
chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']
for chrom in chr_list:
	row_list[chrom]={}
	case_list[chrom]={}
	chrom_size=chrom_size_dic[chrom]
	for i in range(0, chrom_size//bin_size+1):
		row_list[chrom][i]={}
		case_list[chrom][i]={}

in_file=f_list.readline().strip()
while in_file:
	lines=open(in_file)
	sample_id=id_list.readline().strip()
	print('Start making list of '+sample_id)
	in_line=lines.readline().strip()
	while in_line:
		if in_line[0:4]=='#CHR':
			head=in_line+'\tcasNo\tcase_list'
		elif in_line[0]=='#':
			'blank'
		else:
			in_indi=in_line.split('\t')
			idx='\t'.join(in_indi[0:2]+in_indi[3:5])
			chr1=in_indi[0]
			if chr1 not in chr_list:
				in_line=lines.readline().strip()
				continue
			pos1=int(in_indi[1])
			if idx in case_list[chr1][pos1//bin_size].keys():
				prev_case=case_list[chr1][pos1//bin_size][idx]
				case_list[chr1][pos1//bin_size][idx]=prev_case+';'+sample_id
			else:
				case_list[chr1][pos1//bin_size][idx]=sample_id
				row_list[chr1][pos1//bin_size][idx]=in_line
		in_line=lines.readline().strip()
	in_file=f_list.readline().strip()

print('Start sorting and writing')

out_file.write(head+'\n')
for chrom in chr_list:
	print(chrom)
	chrom_size=chrom_size_dic[chrom]
	for i in range(0, chrom_size//bin_size +1):
		sorted_list=[]
		for idx in row_list[chrom][i].keys():
			caseNo=len(case_list[chrom][i][idx].split(';'))
			if caseNo == InCaseNum:
				sorted_list.append([int(idx.split('\t')[1]),caseNo,idx])
		sorted_list.sort(key=itemgetter(0,1))
		for [pos1, caseNo, idx] in sorted_list:
			out_file.write(row_list[chrom][i][idx]+'\t'+str(caseNo)+'\t'+case_list[chrom][i][idx]+'\n')
