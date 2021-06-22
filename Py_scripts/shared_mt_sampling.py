#arg1: input
#arg2: column number of nr3info
#arg3: counted number
#arg4: the number of sampling


#191224 column number change
#200407 edit to get column number

import sys
in_file=open(sys.argv[1])
ncol_nr = int(sys.argv[2])
total_n=int(sys.argv[3])
max_n=int(sys.argv[4])

ncol_nr = ncol_nr -1
out_file=open(sys.argv[1]+'.sam'+str(max_n)+'.txt','w')
#header_list=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','Sample1','ref_readN','var_readN','VAF_pct','nr3ids']
in_line=in_file.readline().strip()
sum_dic={}
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		info='\t'.join(in_indi[0:12])
		refn=int(in_indi[10])
		varn=int(in_indi[11])
		vaf=round(varn*100/float(refn+varn),2)
		nrinfo=in_indi[ncol_nr]
		nr3ids=';'.join(nrinfo.split(';')[4:])
		nr1=int(nrinfo.split(';')[1])
		nr3=int(nrinfo.split(';')[3])
		if nr3==1 or nr1 ==total_n:  # prefilter
			'blank'
		else:
			if nr3ids not in sum_dic.keys():
				sum_dic[nr3ids]=[]
			if len(sum_dic[nr3ids]) < max_n:
				sum_dic[nr3ids].append(in_line)
	in_line=in_file.readline().strip()

n=0
for nr3ids in sum_dic.keys():
	for info in sum_dic[nr3ids]:
		n +=1
		out_file.write(info+'\n')

print(n)
