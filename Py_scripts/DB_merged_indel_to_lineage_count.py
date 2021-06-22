#arg1: input
#arg2: column number of sample_ids (e.g. 7)
#arg3: column number of blood info (e.g. 8)
#arg4: DB_id (e.g. DB2)
#arg5: /path/to/meta_data(e.g. '/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200109.txt')
#arg6: column_number_of lineage_id in meta_data (e.g. 14)
#arg7: /path/to/forcedly_assigned_lineage (e.g.'/home/users/sypark/00_Project/06_LineageTracing/meta_data/forcedly_assigned_lineage.txt') (format: DB_id \t var_id \t lineage_id)


#Using predefined lineage id information and merged indel list, make the count table or full list of point mutation of each lineage
#2019-12-12 header printing error of empty line correction, add blood vaf information
#2020-03-25 modified for indels
#2020-03-29 remove -p option in mkdir perLineage

import sys,os
n_col=int(sys.argv[2])
b_col=int(sys.argv[3])
this_DB=sys.argv[4]
path_meta=sys.argv[5]
l_col=int(sys.argv[6])
path_assign=sys.argv[7]


id_dic={}
NA_sample_list=[]
sample_list=[]
id_file=open(path_meta)
id_line=id_file.readline().strip()
while id_line:
	id_indi=id_line.split('\t')
	DBid=id_indi[0]
	sample_name=id_indi[1]
	l_id=id_indi[l_col-1] # check
	if DBid == this_DB:
		if l_id == 'NA':
			NA_sample_list.append(sample_name)
			id_line=id_file.readline().strip()
			continue
		sample_list.append(sample_name)
		l_id_indi=l_id.split('-')
		for i in range(1,len(l_id_indi)+1):
			this_l_id='-'.join(l_id_indi[0:i])
			if this_l_id not in id_dic.keys():
				id_dic[this_l_id]={}
				id_dic[this_l_id]['samples']=[]
				id_dic[this_l_id]['indels']=[]
				id_dic[this_l_id]['blood_dp']=0
				id_dic[this_l_id]['blood_var']=0
			id_dic[this_l_id]['samples'].append(sample_name)
	id_line=id_file.readline().strip()
id_file.close()

as_file=open(path_assign)
as_line=as_file.readline().strip()
as_dic={}
while as_line:
	as_indi=as_line.split('\t')
	DBid = as_indi[0]
	var_id=as_indi[1]
	lineage_id=as_indi[2]
	if DBid == this_DB:
		as_dic[var_id]=lineage_id
	as_line=as_file.readline().strip()
as_file.close()

in_file=open(sys.argv[1])
mis_file=open(this_DB+'_Possible_missing_lineages.txt','w')
mis_file2=open(this_DB+'_Possible_cooccured.txt','w')
in_line=in_file.readline().strip()
header_line=''
while in_line:
	if in_line[0:4]=='#CHR':
		header_line=in_line
	else:
		in_indi=in_line.split('\t')
		var_id=':'.join([in_indi[0], in_indi[1], in_indi[3], in_indi[4]])
		nr3ids=in_indi[n_col-1].split(';') #check
		nr3ids=list(set(nr3ids) & set(sample_list))
		
		blood_ref=int(in_indi[b_col-1].split(';')[0])
		blood_var=int(in_indi[b_col-1].split(';')[1])
		blood_dp=blood_ref + blood_var
		if var_id in as_dic.keys():
			id_dic[as_dic[var_id]]['indels'].append(in_line)
			id_dic[as_dic[var_id]]['blood_dp'] += blood_dp
			id_dic[as_dic[var_id]]['blood_var'] += blood_var
			in_line=in_file.readline().strip()
			continue
		if len(nr3ids) == 0:
			in_line=in_file.readline().strip()
			continue
		matched='off'
		for l_id in id_dic.keys():
			if set(id_dic[l_id]['samples']) == set(nr3ids):
				id_dic[l_id]['indels'].append(in_line)	
				id_dic[l_id]['blood_dp'] += blood_dp
				id_dic[l_id]['blood_var'] += blood_var
				matched='on'
				break
		if matched=='off':
			print(in_line)
			matching_nums=[]
			matching_ids=[]
			matching_samples=[]
			for l_id in id_dic.keys():
				if set(id_dic[l_id]['samples']) <= set(nr3ids):
					id_dic[l_id]['indels'].append(in_line)
					id_dic[l_id]['blood_dp'] += blood_dp
					id_dic[l_id]['blood_var'] += blood_var
					matching_nums.append(len(id_dic[l_id]['samples']))
					matching_ids.append(l_id)
					matching_samples.append(';'.join(id_dic[l_id]['samples']))
			if sorted(matching_nums, reverse = True)[1] >= 2:
				mis_file.write('##'+in_line+'\n')
				mis_file.write(';'.join(nr3ids)+'\n')
				mis_file.write('\t'.join(matching_ids)+'\n')
				mis_file.write('\t'.join(matching_samples)+'\n\n')
			else:
				mis_file2.write('##'+in_line+'\n')
				mis_file2.write(';'.join(nr3ids)+'\n')
				mis_file2.write('\t'.join(matching_ids)+'\n')
				mis_file2.write('\t'.join(matching_samples)+'\n\n')
				
	in_line=in_file.readline().strip()

	
out_file=open(this_DB+'_Lineage_count_table.txt','w')
out_file.write('##input_file = '+sys.argv[1]+'\n')
ct_header=['#lineage_id','samples','n_samples','n_indels', 'total_blood_dp','total_blood_var','mean_blood_VAFpct']
out_file.write('\t'.join(ct_header)+'\n')
for l_id in id_dic.keys():
	if len(id_dic[l_id]['indels']) == 0:
		continue
#		print('ERROR: no indels was assigned to the lineage. Exiting..')
#		print(l_id)
#		print(id_dic[l_id]['samples'])
#		sys.exit(1)
	blvafpct=id_dic[l_id]['blood_var']*100/float(id_dic[l_id]['blood_dp'])
	info_list=[l_id, ';'.join(id_dic[l_id]['samples']),str(len(id_dic[l_id]['samples'])), str(len(id_dic[l_id]['indels'])), str(id_dic[l_id]['blood_dp']), str(id_dic[l_id]['blood_var']),str(blvafpct)]
	out_file.write('\t'.join(info_list)+'\n')
out_file.close()

os.system('mkdir perLineage')

for l_id in id_dic.keys():
	out_file=open('./perLineage/'+l_id+'.indels.txt','w')
	out_file.write(header_line+'\n')
	out_file.write('\n'.join(id_dic[l_id]['indels'])+'\n')
	out_file.close()

#high_burden_list=['L1-2-1-3-2-1', 'L1-3-2-1-2','L2-3-1-1']
#out_file=open('DB6_shared_in_2more_except_highs_indels.txt','w')
#out_file.write('\t'.join(header_line.split('\t')[0:8])+'\tLineage_id\n')
#for l_id in id_dic.keys():
#	if len(id_dic[l_id]['samples']) >=2 and l_id not in high_burden_list:
#		for in_line in id_dic[l_id]['indels']:
#			in_indi=in_line.split('\t')
#			out_file.write('\t'.join(in_indi[0:8])+'\t'+l_id+'\n')
#out_file.close()
