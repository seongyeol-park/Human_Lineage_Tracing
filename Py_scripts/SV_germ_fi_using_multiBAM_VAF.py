#args1: input file
#args2: column numbers with mVAF info (e.g. 43,44,45,46)

import sys

print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.g_fi','w')
n_cols_list = [int(n)-1 for n in sys.argv[2].split(',')]

n=0; m=0
in_line = in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		n = n +1
		in_indi=in_line.split('\t')
		this_res=[]
		for i in n_cols_list:
			all_dis=int(in_indi[i].split(';')[2])
			this_res.append(all_dis)
		if min(this_res) > 0: #germline
			'blank'
		else:
			m = m+1
			out_file.write(in_line+'\n')
	in_line = in_file.readline().strip()

print(n)
print(m)
