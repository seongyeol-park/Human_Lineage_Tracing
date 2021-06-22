import sys
in_file=open(sys.argv[1])
target_id = sys.argv[2] 
targetNo = int(sys.argv[3]) # include less than this number

out_file=open(sys.argv[1]+'.inc'+target_id+'.less'+str(targetNo),'w')


in_line= in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		caseNo=int(in_indi[10])
		caselist=in_indi[11]
		if caseNo < targetNo and target_id in caselist:
			out_file.write(in_line+'\n')
	in_line = in_file.readline().strip()

