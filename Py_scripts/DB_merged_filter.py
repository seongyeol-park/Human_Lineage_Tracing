#Arg1: input
#Arg2: false snv list
#Arg3: nr1 cutoff
#Arg4: column_number_for_its_numb
#Arg5: column_numbre_for_nr_info


#191015 copied from DB6_merged_filter.py
#210325 space*8 -> tab, make column number as argument

import sys
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.fi','w')

f_file=open(sys.argv[2])
sn1_co=int(sys.argv[3])
col_its = int(sys.argv[4])
col_nr = int(sys.argv[5])

col_its = col_its -1
col_nr = col_nr -1


f_line=f_file.readline().strip()
false_list=[]
while f_line:
	false_list.append(f_line)
	f_line=f_file.readline().strip()


in_line=in_file.readline().strip()
n=m=0
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		n +=1
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=in_indi[1]
		refnt=in_indi[3]
		altnt=in_indi[4]

		its_num=in_indi[col_ids]
		sn_info=in_indi[col_nr]

		idx = '\t'.join([chr1, pos1, refnt, altnt])
		its_num=int(its_num)
		sn1=int(sn_info.split(';')[1])
		sn3=int(sn_info.split(';')[3])
		sn3_list=sn_info.split(';')[4:]

		if sn1-sn3 < 5 and sn3-its_num < 5 and sn1 < sn1_co and idx not in false_list:
			m+=1
			out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()
print(n)
print(m)
