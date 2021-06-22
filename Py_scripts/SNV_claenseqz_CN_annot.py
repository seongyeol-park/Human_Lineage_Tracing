#this script was made 2018-06-14

import sys,os
print(sys.argv[1])
in_file=open(sys.argv[1])  # SNV file
ref_file=open(sys.argv[2])  # sequenza clean_segments file
out_file=open(sys.argv[1]+'.seqzcn','w') 
ref_line=ref_file.readline().strip()
ref_line=ref_file.readline().strip() # pass 1st row
ref_list=[]
while ref_line:
	ref_indi=ref_line.split('\t')
	chr1=ref_indi[0]
	start1=int(ref_indi[1])
	end1=int(ref_indi[2])
	absCN=ref_indi[3]
	Aall=ref_indi[4]
	Ball=ref_indi[5]
	ref_list.append([chr1,start1,end1,absCN,Aall,Ball])
	ref_line=ref_file.readline().strip()

in_line=in_file.readline().strip()
while in_line:
	if in_line[0:6]=='#CHROM':
		header_info=['totCN','majCN','minCN']
		out_file.write(in_line+'\t'+'\t'.join(header_info)+'\n')
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		info1=''
		upst1_dist=1000000000
		downst1_dist=1000000000
		upst1='';downst1=''
		chr1=in_indi[0]
		pos1=int(in_indi[1])
		if chr1=='Y':
			info1='.\t.\t.'
			in_line=in_file.readline().strip()
			continue

		for [ch,st,en,cn,aa,ba] in ref_list:
			if ch==chr1 and en <pos1 and (pos1-en) < upst1_dist:
				upst1='\t'.join([cn,aa,ba])
				upst1_dist=pos1-en
			elif ch==chr1 and pos1>=st and pos1<=en:
				info1='\t'.join([cn,aa,ba])
			elif ch==chr1 and pos1<st and (st-pos1) < downst1_dist:
				downst1='\t'.join([cn,aa,ba])
				downst1_dist=st-pos1
			if info1!='':
				break
				
		if info1=='':
			if upst1_dist < downst1_dist:
				info1=upst1
			elif upst1_dist >= downst1_dist:
				info1=downst1

		out_file.write(in_line+'\t'+info1+'\n')
	in_line=in_file.readline().strip()
