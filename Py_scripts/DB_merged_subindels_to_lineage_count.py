#arg1: input
#arg2: column number of calledid delimited by ; (e.g. 6)
#arg3: column number of blood info (e.g. 7)
#arg4: DB_id (e.g. DB2)
#arg5: /path/to/meta_data(e.g. '/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_20XXXX.txt')
#arg6: column_number_of lineage_id in meta_data (e.g. 15)
#arg7: /path/to/forcedly_assigned_lineage (e.g.'/home/users/sypark/00_Project/06_LineageTracing/meta_data/forcedly_assigned_lineage_cindel.txt') (format: DB_id \t var_id \t lineage_id)

#Using predefined lineage id information and merged substitution list, make the count table or full list of point mutation of each lineage
#2019-12-12 header printing error of empty line correction, add blood vaf information
#2020-03-25 write possible co-occur list; algorithm for finding possible was updated
#2020-04-03 update algorithm for listing missed lineage: deduplication of ids, add RM in forcedly_assigned_list; remove -p option in mkdir
#2020-04-09 add n_subs and n_indels separately
#2020-05-22 only include samples which are marked as "included" in meta data.
#2020-06-19 multiple lineage id can be entered in forcedly assigned vars

import sys,subprocess
n_col=int(sys.argv[2])
b_col=int(sys.argv[3])
this_DB=sys.argv[4]
path_meta=sys.argv[5]
l_col=int(sys.argv[6])
path_assign=sys.argv[7]

res = subprocess.Popen('mkdir perLineage', shell = True)
if res.wait() != 0:
	print('Error. Exit')
	sys.exit(1)

id_dic={}
NA_sample_list=[]
sample_list=[]
id_file=open(path_meta)
id_line=id_file.readline().strip()
while id_line:
	id_indi=id_line.split('\t')
	DBid=id_indi[0]
	sample_name=id_indi[1]
	included=id_indi[2]
	if included != 'Y':
		id_line = id_file.readline().strip()
		continue
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
				id_dic[this_l_id]['pointmts']=[]
				id_dic[this_l_id]['subs']=[]
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
		refnt = in_indi[3]
		altnt = in_indi[4]
		if len(refnt) == len(altnt):
			mttype = 'subs'
		else:
			mttype = 'indels'
		nr3ids=in_indi[n_col-1].split(';')  #check
		nr3ids=list(set(nr3ids) & set(sample_list))
			

		blood_dp=int(in_indi[b_col-1].split(';')[0])
		blood_var=int(in_indi[b_col-1].split(';')[1])
		if var_id in as_dic.keys():
			for lid in as_dic[var_id].split(','):
				id_dic[lid]['pointmts'].append(in_line)
				if mttype == 'subs':
					id_dic[lid]['subs'].append(in_line)
				else:
					id_dic[lid]['indels'].append(in_line)
				id_dic[lid]['blood_dp'] += blood_dp
				id_dic[lid]['blood_var'] += blood_var
			in_line=in_file.readline().strip()
			continue

		if len(nr3ids) == 0:
			in_line=in_file.readline().strip()
			continue
		matched='off'
		for l_id in id_dic.keys():
			if set(id_dic[l_id]['samples']) == set(nr3ids):
				id_dic[l_id]['pointmts'].append(in_line)	
				if mttype == 'subs':
					id_dic[l_id]['subs'].append(in_line)
				else:
					id_dic[l_id]['indels'].append(in_line)
				id_dic[l_id]['blood_dp'] += blood_dp
				id_dic[l_id]['blood_var'] += blood_var
				matched='on'
				break
		if matched=='off':
			matching_nums=[]
			matching_ids=[]
			matching_id_levels=[]
			matching_samples=[]
			matching_n = 0
			for l_id in id_dic.keys():
				if set(id_dic[l_id]['samples']) <= set(nr3ids):
					id_dic[l_id]['pointmts'].append(in_line)
					if mttype == 'subs':
						id_dic[l_id]['subs'].append(in_line)
					else:
						id_dic[l_id]['indels'].append(in_line)
					id_dic[l_id]['blood_dp'] += blood_dp
					id_dic[l_id]['blood_var'] += blood_var
					matching_nums.append(len(id_dic[l_id]['samples']))
					matching_ids.append(l_id)
					matching_id_levels.append(len(l_id.split('-')))

			dedup_matching_ids=[]
			for l_id1 in matching_ids:
				dup='off'
				for l_id2 in matching_ids:
					if len(l_id2) >= len(l_id1):
						continue
					if l_id1[0:len(l_id2)] == l_id2:
						dup ='on'
				if dup == 'off':
					dedup_matching_ids.append(l_id1)
					matching_samples.append(';'.join(id_dic[l_id1]['samples']))
					matching_n = matching_n + len(id_dic[l_id1]['samples'])
			core_id_list=[]
			for l_id in dedup_matching_ids:
				core_id = '-'.join(l_id.split('-')[0:min(matching_id_levels)-1])
				core_id_list.append(core_id)

			if len(list(set(core_id_list))) == 1:
				mis_file.write('##'+in_line+'\n')
				mis_file.write(';'.join(nr3ids)+'\n')
				mis_file.write('\t'.join(dedup_matching_ids)+'\n')
				mis_file.write('\t'.join(matching_samples)+'\n')
				mis_file.write('No. of matched samples : '+str(matching_n)+'\n\n')
			else:
				mis_file2.write('##'+in_line+'\n')
				mis_file2.write(';'.join(nr3ids)+'\n')
				mis_file2.write('\t'.join(dedup_matching_ids)+'\n')
				mis_file2.write('\t'.join(matching_samples)+'\n')
				mis_file2.write('No. of matched samples : '+str(matching_n)+'\n\n')
				
#			sorted(matching_nums, reverse = True)[1]
#			if sorted(matching_nums, reverse = True)[1] >= 2:
#				mis_file.write('##'+in_line+'\n')
#				mis_file.write(';'.join(nr3ids)+'\n')
#				mis_file.write('\t'.join(matching_ids)+'\n')
#				mis_file.write('\t'.join(matching_samples)+'\n')
#			else:
#				mis_file2.write('##'+in_line+'\n')
#				mis_file2.write(';'.join(nr3ids)+'\n')
#				mis_file2.write('\t'.join(matching_ids)+'\n')
#				mis_file2.write('\t'.join(matching_samples)+'\n')
	in_line=in_file.readline().strip()

	
out_file=open(this_DB+'_Lineage_count_table.txt','w')
out_file.write('##input_file = '+sys.argv[1]+'\n')
ct_header=['#lineage_id','samples','n_samples','n_pointmt','n_subs','n_indels','total_blood_dp','total_blood_var','mean_blood_VAFpct']
out_file.write('\t'.join(ct_header)+'\n')
for l_id in id_dic.keys():
	if len(id_dic[l_id]['pointmts']) == 0:
		print('ERROR: no pointmts was assigned to the lineage. Exiting..')
		print(l_id)
		print(id_dic[l_id]['samples'])
		sys.exit(1)
	blvafpct=id_dic[l_id]['blood_var']*100/float(id_dic[l_id]['blood_dp'])
	info_list=[l_id, ';'.join(id_dic[l_id]['samples']),str(len(id_dic[l_id]['samples'])), str(len(id_dic[l_id]['pointmts'])), str(len(id_dic[l_id]['subs'])), str(len(id_dic[l_id]['indels'])), str(id_dic[l_id]['blood_dp']), str(id_dic[l_id]['blood_var']),str(blvafpct)]
	out_file.write('\t'.join(info_list)+'\n')
out_file.close()

for l_id in id_dic.keys():
	out_file=open('./perLineage/'+l_id+'.pointmts.txt','w')
	out_file.write(header_line+'\n')
	out_file.write('\n'.join(id_dic[l_id]['pointmts'])+'\n')
	out_file.close()

