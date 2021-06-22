#args1: id \t input file list
#args2: output file name

import sys,os
fn_file=open(sys.argv[1])  # id \t file (chr1 \t pos1 \t chr2 \t pos2 \t mh \t ori... )
fn_line=fn_file.readline().strip()
out_file=open(sys.argv[2],'w')
header_list=[]
line_list=[]
while fn_line:
	print(fn_line)
	fn_indi=fn_line.split('\t')
	in_file=open(fn_indi[1])
	sampleid=fn_indi[0]
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0:4]=='#CHR':
			header_list.append(in_line)
		elif in_line[0]=='#':
			'blank'
		else:	
			line_list.append(in_line+'\t'+sampleid)
		in_line=in_file.readline().strip()
	fn_line=fn_file.readline().strip()

header_list=list(set(header_list))
if len(header_list) > 1:
	print('ERROR:Multiple types of Header.exit')
	print(header_list)
	sys.exit(1)


out_file.write(header_list[0]+'\tsampleid\n')
for line in line_list:
	out_file.write(line+'\n')
out_file.close()

os.system('python /home/users/sypark/01_Python_files/strucural_variation/SV_sorting.py '+sys.argv[2])
os.system('rm '+sys.argv[2])
			
