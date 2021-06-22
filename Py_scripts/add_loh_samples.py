#arg1: input
#arg2: id \t sequenza segments.txt.clean
#arg3: gender (M or F)

#200106: gender setting added
#200907: error correction in gender == 'M': add chr1 not in ['X','Y'] and

import sys

seg_fn=open(sys.argv[2])
gender=sys.argv[3]
seg_line=seg_fn.readline().strip()
ref_dic={}
while seg_line:
	seg_indi=seg_line.split('\t')
	print(seg_line)
	sample_id=seg_indi[0]
	if sample_id not in ref_dic.keys():
		ref_dic[sample_id]={}
	ref_file=open(seg_indi[1])
	ref_line=ref_file.readline().strip()
	while ref_line:
		if ref_line[0]=='#':
			'blank'
		else:
			ref_indi=ref_line.split('\t')
			chr1=ref_indi[0]
			start_pos=ref_indi[1]
			end_pos=ref_indi[2]
			tcn=ref_indi[3]
			mcn=ref_indi[5]
			if chr1 not in ref_dic[sample_id].keys():
				ref_dic[sample_id][chr1]={}
			ref_dic[sample_id][chr1][start_pos+'\t'+end_pos]=[int(tcn),int(mcn)]
		ref_line=ref_file.readline().strip()
	seg_line=seg_fn.readline().strip()

in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.loh_id','w')
in_line=in_file.readline().strip()
while in_line:
	if in_line[0:4]=='#CHR':
		out_file.write(in_line+'\tLOH_ids\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		if gender == 'F' and chr1 == 'Y':
			in_line=in_file.readline().strip()
			continue
		pos1=int(in_indi[1])
		loh_list=[]
		for sample in ref_dic.keys():
			for info in ref_dic[sample][chr1].keys():
				start_pos=int(info.split('\t')[0])
				end_pos=int(info.split('\t')[1])
				mcn=ref_dic[sample][chr1][info][1]
				if gender == 'F':
					if pos1 >= start_pos and pos1 <= end_pos and mcn == 0: 
						loh_list.append(sample)
				elif gender == 'M':
					if chr1 in ['X','Y'] and pos1 >= start_pos and pos1 <= end_pos and tcn == 0:
						loh_list.append(sample)
					elif chr1 not in ['X','Y'] and pos1 >= start_pos and pos1 <= end_pos and mcn == 0:
						loh_list.append(sample)
		if len(loh_list)==0:
			out_file.write(in_line+'\tNA\n')
		else:
			out_file.write(in_line+'\t'+';'.join(loh_list)+'\n')
	in_line=in_file.readline().strip()
