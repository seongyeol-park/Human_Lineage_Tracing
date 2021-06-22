#2021-03-13 make this compatible with python3

import sys, gzip
if sys.argv[1][-3:]=='.gz':
	in_file=gzip.open(sys.argv[1],'rt')
else:
	in_file=open(sys.argv[1]) #non-redundant sorted vcf1

if sys.argv[2][-3:]=='.gz':
	ref_file=gzip.open(sys.argv[2],'rt')
else:
	ref_file=open(sys.argv[2]) #non redundant sorted vcf2
sampleid= sys.argv[3]
caller1=sys.argv[4]
caller2=sys.argv[5]
print(sampleid)

out_file3=open(sys.argv[3]+'.'+caller1+'_'+caller2+'_common','w')
in_line=in_file.readline().strip()
while in_line[0]=='#':
	if in_line[0:6]=='#CHROM':
		header1=in_line
	in_line=in_file.readline().strip()

ref_line=ref_file.readline().strip()
while ref_line[0]=='#':
	if ref_line[0:6]=='#CHROM':
		header2=ref_line
	ref_line=ref_file.readline().strip()


#out_file1.write(header1+'\n')
#out_file2.write(header2+'\n')
out_file3.write(header1+'\n')

common=0
in_only=0
ref_only=0
in_dis=0
ref_dis=0

discard_list=['GL','hs','NC']
prv_chr1=0; prv_chr2=0
while in_line or ref_line:	
	if in_line[0:2] in discard_list:
		in_dis=in_dis+1
		in_line=in_file.readline().strip()
	elif ref_line[0:2] in discard_list:
		ref_dis=ref_dis+1
		ref_line=ref_file.readline().strip()
	elif ref_line=='':
		in_only=in_only+1
		#out_file1.write(in_line+'\n')
		in_line=in_file.readline().strip()
	elif in_line=='':
		ref_only=ref_only+1
		#out_file2.write(ref_line+'\n')
		ref_line=ref_file.readline().strip()
	else:
		in_indi=in_line.split('\t')
		ref_indi=ref_line.split('\t')
		in_idx=in_indi[0]+'\t'+in_indi[1]+'\t'+in_indi[3]+'\t'+in_indi[4]
		ref_idx=ref_indi[0]+'\t'+ref_indi[1]+'\t'+ref_indi[3]+'\t'+ref_indi[4]
		chr1=int(((in_indi[0].replace('X','23')).replace('Y','24')).replace('MT','25'))
		pos1=int(in_indi[1])
		chr2=int(((ref_indi[0].replace('X','23')).replace('Y','24')).replace('MT','25'))
		pos2=int(ref_indi[1])
		if chr1==chr2:
			if pos1==pos2:
				if in_idx==ref_idx:
					common=common+1
					out_file3.write(in_line+'\n')
					in_line=in_file.readline().strip()
					ref_line=ref_file.readline().strip()
				else:
					in_only=in_only+1
					ref_only=ref_only+1
					#out_file1.write(in_line+'\n')
					in_line=in_file.readline().strip()
					#out_file2.write(ref_line+'\n')
					ref_line=ref_file.readline().strip()
			elif pos1>pos2:
				ref_only=ref_only+1
				#out_file2.write(ref_line+'\n')
				ref_line=ref_file.readline().strip()
			elif pos1<pos2:
				in_only=in_only+1
				#out_file1.write(in_line+'\n')
				in_line=in_file.readline().strip()
		elif chr1>chr2:
			ref_only=ref_only+1
			#out_file2.write(ref_line+'\n')
			ref_line=ref_file.readline().strip()
		elif chr1<chr2:
			in_only=in_only+1
			#out_file1.write(in_line+'\n')
			in_line=in_file.readline().strip()
		if prv_chr1 > chr1:
			print('Fatal error, exiting: Wrong sorting: '+sys.argv[1])
			sys.exit()
		if prv_chr2 > chr2:
			print('Fatal error, exiting, Wrong sorting: '+sys.argv[2])
			sys.exit()
		prv_chr1 = chr1
		prv_chr2 = chr2
	
print('Common= '+str(common))
print('First file only= '+str(in_only))
print('Second file only= '+str(ref_only))
