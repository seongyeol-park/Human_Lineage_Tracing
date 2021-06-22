#arg1: input (id \t path )
#arg2: meta_data_path (sampleid at column2, lineageid at column15)
#arg3: deadbody id (e.g. DB6)

import sys
from operator import itemgetter

fn_file=open(sys.argv[1]) # id \t file_path
meta_file=open(sys.argv[2])
thisDB = sys.argv[3]

var_case_dic={}
list_for_sort=[]
list_for_compare=[]
fn_line = fn_file.readline().strip()
while fn_line:
	print(fn_line)
	fn_indi=fn_line.split('\t')
	sampleid = fn_indi[0]
	filepath = fn_indi[1]
	in_file=open(filepath)
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0:4] == '#CHR':
			header = in_line
		elif in_line[0]=='#':
			'blank'
		else:
			in_indi=in_line.split('\t')
			chr1 = in_indi[0]
			chr1n = int(chr1.replace('X','23').replace('Y','24'))
			pos1 = in_indi[1]
			chr2 = in_indi[2]
			pos2 = in_indi[3]
			mh = in_indi[4]
			ter = in_indi[5]
			svtype = in_indi[6]
			varid = ':'.join([chr1, pos1, chr2, pos2,mh, ter, svtype])
			list_for_compare.append([chr1, int(pos1), chr2, int(pos2.replace('.','0')), ter[0], ter[-1], varid, sampleid])
			if varid in var_case_dic.keys():
				var_case_dic[varid] = var_case_dic[varid]+';'+sampleid
			else:
				list_for_sort.append([chr1n, int(pos1), varid])
				var_case_dic[varid] = sampleid
		in_line = in_file.readline().strip()
	fn_line=fn_file.readline().strip()

out_file1=open('Possible_shared_SVs.txt','w')
n=0
for [chr1, pos1, chr2, pos2, ter1, ter2, varid1, sample1] in list_for_compare:
	for [chr3, pos3, chr4, pos4, ter3, ter4, varid2, sample2] in list_for_compare:
		if varid1 == varid2:
			continue
		else:
			if chr1 == chr3 and chr2 == chr4 and ter1 == ter3 and ter2 == ter4 and abs(pos1-pos3) < 1000 and abs(pos2-pos4) < 1000:
				n  = n+1
				out_file1.write('#'+str(n)+'\n')
				out_file1.write('\t'.join(varid1.split(':'))+'\t'+sample1+'\n')
				out_file1.write('\t'.join(varid2.split(':'))+'\t'+sample2+'\n')

meta_line = meta_file.readline().strip()
meta_line = meta_file.readline().strip() # pass 1st row
SLdic={}
while meta_line:
	meta_indi = meta_line.split('\t')
	#column number check
	DB_id = meta_indi[0]
	s_id = meta_indi[1]  
	l_id = meta_indi[14]
	if DB_id == thisDB:
		SLdic[s_id]=l_id
	meta_line=meta_file.readline().strip()

out_file2 = open(thisDB+'_Merged_SV_list.txt','w')
header_list=['#CHR1','POS1','CHR2','POS2','mh','terinfo','svtype','sample_id','lineage_id']
out_file2.write('\t'.join(header_list)+'\n')
list_for_sort.sort(key = itemgetter(0,1))
for [chr1n, pos1, varid1] in list_for_sort:
	samples = var_case_dic[varid1]
	n_sample = len(samples.split(';'))
	if n_sample > 1:
		sample_list = samples.split(';')
		ln_list=[]
		ln_len_list=[]
		for sn in sample_list:
			ln_list.append(SLdic[sn])
			ln_len_list.append(len(SLdic[sn]))
		cor_ln_list=[]
		min_len = min(ln_len_list)
		print(ln_list)
		for ln in ln_list:
			cor_ln_list.append(ln[0:min_len-2])
		print(cor_ln_list)
		cor_ln_list = list(set(cor_ln_list))
		if len(cor_ln_list) > 1:
			print('###FATAL ERROR. fail to assign lineage id')
			print('\t'.join(varid1.split(';')))
			print('\t'.join(sample_list))
			sys.exit(1)
		final_lid = cor_ln_list[0]
	else:
		final_lid = SLdic[samples]
	out_file2.write('\t'.join(varid1.split(':'))+'\t'+samples+'\t'+final_lid+'\n')

