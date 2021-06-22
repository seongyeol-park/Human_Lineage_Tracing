#arg1: input
#arg2: false_SV_list

import sys
f_file=open(sys.argv[2])
f_line=f_file.readline().strip()
false_list=[]
while f_line:
	false_list.append(f_line)
	f_line = f_file.readline().strip()


n=0; m=0
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.rmFP','w')
in_line = in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		n = n+1
		in_indi=in_line.split('\t')
		sv_id='_'.join(in_indi[0:4]+[in_indi[5]])
		if sv_id not in false_list:
			m = m+1
			out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()

print(n)
print(m)
