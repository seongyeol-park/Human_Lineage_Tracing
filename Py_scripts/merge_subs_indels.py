#args1: merged substitution list
#args2: merged indel list
#args3: DB_id (e.g. DB2)

#200402: written
#200404 if DB6: column number 

import sys
from operator import itemgetter

sub_file=open(sys.argv[1])
indel_file=open(sys.argv[2])
out_list=[]

sub_line=sub_file.readline().strip()
while sub_line:
	if sub_line[0] == '#':
		sub_line = sub_file.readline().strip()
		continue
	sub_indi = sub_line.split('\t')
	chr1 = sub_indi[0]
	chr1n = int(chr1.replace('X','23').replace('Y','24'))
	pos1 = int(sub_indi[1])
	refnt = sub_indi[3]
	altnt = sub_indi[4]
	if sys.argv[3] == 'DB6':
		nr3info = sub_indi[45]
		bloodinfo = sub_indi[46]
	else:
		nr3info = sub_indi[47]
		bloodinfo = sub_indi[48]

	called_ids = ';'.join(nr3info.split(';')[4:])
	bloodinfo_re = bloodinfo.replace('%','')

	merged_info='\t'.join([chr1, str(pos1), '.', refnt, altnt, called_ids, bloodinfo_re])
	out_list.append([chr1n, pos1, merged_info])
	sub_line = sub_file.readline().strip()
sub_file.close()


indel_line = indel_file.readline().strip()
while indel_line:
	if indel_line[0] == '#':
		indel_line = indel_file.readline().strip()
		continue
	indel_indi = indel_line.split('\t')
	chr1 = indel_indi[0]
	chr1n = int(chr1.replace('X','23').replace('Y','24'))
	pos1 = int(indel_indi[1])
	refnt = indel_indi[3]
	altnt = indel_indi[4]
	called_ids = indel_indi[5]
	bloodinfo = indel_indi[7]

	blref = int(bloodinfo.split(';')[0])
	blvar = int(bloodinfo.split(';')[1])
	blvaf = bloodinfo.split(';')[3]
	bldp = blref + blvar
	bloodinfo_re = ';'.join([str(bldp), str(blvar), blvaf])

	merged_info='\t'.join([chr1, str(pos1), '.', refnt, altnt, called_ids, bloodinfo_re])
	out_list.append([chr1n, pos1, merged_info])
	indel_line = indel_file.readline().strip()
indel_file.close()

out_list.sort(key = itemgetter(0,1))
out_file=open(sys.argv[3] +'_merged_subs_indels.txt','w')
out_file.write('\t'.join(['#CHROM','POS','ID','REF','ALT','sample_ids','blood_dp;blood_var;blood_vafpct'])+'\n')
for [chr1n, pos1, info] in out_list:
	out_file.write(info+'\n')

