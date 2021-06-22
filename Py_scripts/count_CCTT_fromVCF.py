#arg1: input file list
import sys
fn_file=open(sys.argv[1]) 
out_file=open('cctt_count.txt','w')
out_file.write('file_name\tCCTTcount\n')
fn_line=fn_file.readline().strip()
while fn_line:
	print(fn_line)
	prv_chr=''
	prv_pos=0
	prv_ref=''
	prv_alt=''
	cctt=0
	in_file=open(fn_line)
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':
			'blank'
		else:
			in_indi=in_line.split('\t')
			chr1=in_indi[0]
			pos1=int(in_indi[1])
			refnt=in_indi[3][0]
			altnt=in_indi[4][0]
			if prv_chr == chr1 and pos1-prv_pos == 1 and refnt+altnt in ['CT','GA'] and prv_ref+prv_alt in ['CT','GA']:
				cctt +=1
			prv_chr=chr1; prv_pos=pos1; prv_ref=refnt; prv_alt=altnt
		in_line=in_file.readline().strip()
	out_file.write(fn_line+'\t'+str(cctt)+'\n')
	fn_line=fn_file.readline().strip()
