#arg1: input
#arg2: column number of nr info
#arg3: column number of loh id info
#arg4: meta_dt path 
#arg5: deadbody id (e.g. DB2)

#2020-04-07: edited to get column number by arguments; edtied to get meta data path

import sys
print(sys.argv[1])

in_file=open(sys.argv[1])
ncol_nr = int(sys.argv[2])
ncol_loh = int(sys.argv[3])
meta_file = open(sys.argv[4])
thisDB = sys.argv[5]

ncol_nr = ncol_nr -1
ncol_loh = ncol_loh -1


id_list=[]
meta_line = meta_file.readline().strip() # pass 1st row
meta_line = meta_file.readline().strip()
while meta_line:
	meta_indi=meta_line.split('\t')
	dbid = meta_indi[0]
	sampleid = meta_indi[1]
	included = meta_indi[2]
	if dbid == thisDB and included == 'Y':
		id_list.append(sampleid)
	meta_line = meta_file.readline().strip()


out_file=open(sys.argv[1]+'.g_check','w')
in_line=in_file.readline().strip()
n=0; m=0
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\tgermline\n')
	else:
		n +=1
		in_indi=in_line.split('\t')
		nr_info=in_indi[ncol_nr]
		loh_info=in_indi[ncol_loh]
		nr3=int(nr_info.split(';')[3])
		nr3_ids=nr_info.split(';')[4:]
		loh_ids=loh_info.split(';')
		if nr3_ids[0]=='All':
			out_file.write(in_line+'\tY\n')
		elif len(set(id_list) -set(nr3_ids) - set(loh_ids)) == 0:
			out_file.write(in_line+'\tY\n')
		else:
			m +=1
			out_file.write(in_line+'\tN\n')
	in_line = in_file.readline().strip()

print(n)
print(m)
