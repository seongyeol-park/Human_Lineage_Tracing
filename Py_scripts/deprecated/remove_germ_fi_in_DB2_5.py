import sys

OLcols = [46, 47, 48]  # column number of other lineage rasm info
OLcols = [x-1 for x in OLcols]

in_file=open(sys.argv[1])
print(sys.argv[1])
out_file=open(sys.argv[1]+'.g_fi','w')
in_line=in_file.readline().strip()
n=0; m=0
while in_line:
	if in_line[0]== '#':
		out_file.write(in_line+'\n')
	else:
		n = n+1
		in_indi=in_line.split('\t')
		info_list=[]
		for i in OLcols:
			ncons=int(in_indi[i].split(';')[8])
			info_list.append(ncons)
		if min(info_list) > 0:
			'blank'
		else:
			m = m+1
			out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()

print(n)
print(m)
