#arg1: path of input files
#arg2: common variants list

#200309 if chr1 not in chr_list -> continue

import sys, os
from operator import itemgetter

f_list=open(sys.argv[1])  

chrom_size_dic={}
chr_file=open('/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fa.fai')
chr_line=chr_file.readline().strip()
while chr_line:
	chr_indi=chr_line.split('\t')
	chrom_size_dic[chr_indi[0]]=int(chr_indi[1])
	chr_line=chr_file.readline().strip()
var_list={}
bin_size=100000
chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']
for chrom in chr_list:
	var_list[chrom]={}
	chrom_size=chrom_size_dic[chrom]
	for i in range(0, chrom_size//bin_size+1):
		var_list[chrom][i]=[]


var_file=open(sys.argv[2]) 
var_line=var_file.readline().strip()
while var_line:
	if var_line[0]=='#':
		'blank'
	else:
		var_indi=var_line.split('\t')
		chr1=var_indi[0]
		pos1=int(var_indi[1])
		var_list[chr1][pos1//bin_size].append(pos1)
	var_line=var_file.readline().strip()

fn_line=f_list.readline().strip()
while fn_line:
	print(fn_line)
	n=0; m=0
	in_file=open(fn_line)
	out_file=open(fn_line+'.cmn_fi','w')
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':
			out_file.write(in_line+'\n')
		else:
			n +=1
			in_indi=in_line.split('\t')
			chr1=in_indi[0]
			if chr1 not in chr_list:
				in_line=in_file.readline().strip()
				continue
			pos1=int(in_indi[1])
			if pos1 in var_list[chr1][pos1//bin_size]:
				'blank'
			else:
				m +=1
				out_file.write(in_line+'\n')
		in_line=in_file.readline().strip()
	print(n)
	print(m)
	fn_line=f_list.readline().strip()
