#Arg1 filnal file list with annovar and nr3ids
#Arg2 output name

import sys,collections
fn_file=open(sys.argv[1])  # final files
out_file=open(sys.argv[2],'w')

col_loca=24
col_exonic=27
col_nr3=50

col_loca-=1;col_exonic-=1;col_nr3-=1

fn_line=fn_file.readline().strip()
ids_dic={}
while fn_line:
	print(fn_line)
	in_file=open(fn_line)
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':	
			'blank'
		else:
			in_indi=in_line.split('\t')
			ids=';'.join(in_indi[col_nr3].split(';')[4:])
			if ids not in ids_dic.keys():
				ids_dic[ids]={}
				ids_dic[ids]["loca"]=[]
				ids_dic[ids]["exonic"]=[]
			loca=in_indi[col_loca].split(';')[0]   
			exonic=in_indi[col_exonic]
			ids_dic[ids]["loca"].append(loca)
			if loca=='exonic':
				ids_dic[ids]["exonic"].append(exonic)
		in_line=in_file.readline().strip()
	fn_line=fn_file.readline().strip()


loca_classes=["exonic","splicing","intronic","UTR3","UTR5","ncRNA_exonic","ncRNA_splicing","ncRNA_intronic","upstream","downstream","intergenic"]
exonic_classes=["nonsynonymous SNV","stopgain","stoploss","unknown","synonymous SNV"]
header_list=['#sample','n_sample','n_subs']+loca_classes+exonic_classes
out_file.write('\t'.join(header_list)+'\n')
for ids in ids_dic.keys():
	n_sample=len(ids.split(';'))
	num_mts=len(ids_dic[ids]["loca"])
	loca_dic=collections.Counter(ids_dic[ids]["loca"])
	exonic_dic=collections.Counter(ids_dic[ids]["exonic"])
	out_file.write(ids+'\t'+str(n_sample)+'\t'+str(num_mts))
	for c in loca_classes:
		if c not in loca_dic.keys():
			out_file.write('\t0')
		else:
			out_file.write('\t'+str(loca_dic[c]))
	for c in exonic_classes:
		if c not in exonic_dic.keys():
			out_file.write('\t0')
		else:
			out_file.write('\t'+str(exonic_dic[c]))
	out_file.write('\n')
	
	
	
	
	
