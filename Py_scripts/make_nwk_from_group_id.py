#Arg1: input file with group id and number of mutations


import sys
in_file=open(sys.argv[1]) # count table
col_nsub=4
col_gid=1

col_nsub-=1
col_gid-=1
out_file=open(sys.argv[1]+'.nwk','w')
in_line=in_file.readline().strip()
nsub_dic={}
g_list=[]
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		nsub=int(in_indi[col_nsub])
		gid=in_indi[col_gid]
		g_list.append(gid)
		nsub_dic[gid]=nsub
	in_line=in_file.readline().strip()

print(len(nsub_dic))

branch_num=0
for gid in g_list:
	if len(gid.split('-')) > branch_num:
		branch_num=len(gid.split('-'))
	
nwk=''
#ancestor registration
cur_list=[]
for gid in g_list:
	if len(gid.split('-'))==1:
		cur_list.append(gid)
cur_list.sort()
nwk='('+','.join(cur_list)+')'
for gid in cur_list:
	nwk=nwk.replace(gid, gid+':'+str(nsub_dic[gid]))

#daughter registration
for n in range(2,branch_num+1):
	print(n)
	print(nwk)
	cur_list=[]
	cur_dic={}
	for gid in g_list:
		if len(gid.split('-'))==n:
			cur_list.append(gid)
	cur_list.sort()
	for gid in cur_list:
		mother='-'.join(gid.split('-')[0:-1])
		if mother not in cur_dic.keys():
			cur_dic[mother]=[]
		cur_dic[mother].append(gid)
	for mother in cur_dic.keys():
		nwk=nwk.replace(mother, '('+','.join(cur_dic[mother])+')'+mother)
	for gid in cur_list:
		nwk=nwk.replace(gid, gid+':'+str(nsub_dic[gid]))

out_file.write(nwk+'L0;\n')
