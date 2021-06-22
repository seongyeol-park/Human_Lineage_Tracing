#Arg1: input
#Arg2: correction list (e.g. '/home/users/sypark/00_Project/06_LineageTracing/meta_data/variant_manual_correction.txt')
#Arg3: deadbody id (e.g. DB2)


import sys

print(sys.argv[1])
in_file=open(sys.argv[1])
cs_file=open(sys.argv[2])
thisDB = sys.argv[3]

out_file=open(sys.argv[1]+'.edit','w')

filter_list=[]
cs_line = cs_file.readline().strip()
while cs_line:
	cs_indi = cs_line.split('\t')
	dbid = cs_indi[0]
	varid = cs_indi[1]
	action = cs_indi[2]
	if dbid == thisDB and action == 'rm':
		filter_list.append(varid)
	cs_line = cs_file.readline().strip()

n=0; m=0
in_line = in_file.readline().strip()
while in_line:
	if in_line[0] == '#':
		out_file.write(in_line + '\n')
	else:
		n = n +1
		in_indi = in_line.split('\t')
		chr1 = in_indi[0]
		pos1 = in_indi[1]
		refnt = in_indi[3]
		altnt = in_indi[4]
		
		var_id = ':'.join([chr1, pos1, refnt, altnt])
		
		if var_id not in filter_list:
			m = m +1
			out_file.write(in_line+'\n')
	in_line = in_file.readline().strip()

print(n)
print(m)


