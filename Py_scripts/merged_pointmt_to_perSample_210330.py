#Arg1: input
#Arg2: column number of sample ids (id1;id2;...) (e.g. 6)
#Arg3: DB id (e.g. DB2)
#Arg4: lineage table (e.g. 'DB2_Lineage_count_table.txt') (format: Lineage_id \t sample_ids \t ...)
#Arg5: forcedly assigned lineage (e.g. '/home/users/sypark/00_Project/06_LineageTracing/meta_data/forcedly_assigned_lineage.txt')
#Arg6: type of variants (snv , indel, pointmt)

#191213 generalization
#200109 add arguments and use forcedly assigned lineage
#200329 remove -p option in mkdir perSample
#200403 "rm" processing in forcedly_assigned_list
#200619 multiple lid can be entered in forcedly assinged lineage
#210330 error_correction: if variant is exist in forcedly assigned lineage, sample list in vcf will be ignored.

import sys,subprocess
in_file=open(sys.argv[1])
col_sids = int(sys.argv[2])-1
this_DB = sys.argv[3]
path_li = sys.argv[4]
path_assign = sys.argv[5]
var_type = sys.argv[6]


li_file = open(path_li)
li_line = li_file.readline().strip()
ls_dic={}
while li_line:
	if li_line[0]=='#':
		'blank'
	else:
		li_indi = li_line.split('\t')
		l_id = li_indi[0]
		s_id = li_indi[1]
		ls_dic[l_id]=s_id.split(';')
	li_line=li_file.readline().strip()

as_dic={}
as_file=open(path_assign)
as_line=as_file.readline().strip()
while as_line:
	as_indi=as_line.split('\t')
	db_id=as_indi[0]
	var_id=as_indi[1]
	l_id=as_indi[2]
	if db_id == this_DB and l_id != 'rm':
		as_dic[var_id]=[]
		for indi_lid in l_id.split(','):
			as_dic[var_id]=as_dic[var_id]+ls_dic[indi_lid]
	as_line=as_file.readline().strip()

res_dic={}
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		info='\t'.join(in_indi[0:5])
		var_id=':'.join([in_indi[0], in_indi[1], in_indi[3], in_indi[4]])
		if var_id in as_dic.keys():
			if as_dic[var_id] == 'rm':
				in_line=in_file.readline().strip()
				continue
			sample_list = list(set(as_dic[var_id]))
		else:
			if var_type == 'snv':
				sample_list=in_indi[col_sids].split(';')[4:]
			else:
				sample_list=in_indi[col_sids].split(';')
		for sample in sample_list:
			if sample not in res_dic.keys():
				res_dic[sample]=[]
			res_dic[sample].append(info)
	in_line=in_file.readline().strip()

res = subprocess.Popen('mkdir perSample', shell = True)
if res.wait() != 0:
	print('Error. exit')
	sys.exit(1)


for sample in res_dic.keys():
	print(sample)
	out_file=open('./perSample/'+sample+'.'+var_type,'w')
	header_list=['#CHROM','POS','ID','REF','ALT']
	out_file.write('\t'.join(header_list)+'\n')
	for info in res_dic[sample]:
		out_file.write(info+'\n')
	out_file.close()
	
