#arg1: input
#arg2: chr1 sampleN capculated MPU path  


#191105 copied from DB6_AddSampleNinfo.py, change sampleN and path to argument
#200406 gzip open ('rt') add; add [0] in refnt and altnt for MNVs


import sys,gzip
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.sampleN','w')
in_line=in_file.readline().strip()
prv_chr='0'
nt_list=["A","G","C","T"]
while in_line:
	if in_line[0:4]=='#CHR':
		out_file.write(in_line+'\tdp;Nr1;Nr2;Nr3;Nr3ids\tblood_dp;read;vaf\n')
		in_line=in_file.readline().strip()
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=int(in_indi[1])
		refnt=in_indi[3][0]
		altnt=in_indi[4][0]
		altidx=nt_list.index(altnt)
	
		if chr1 != prv_chr:
			ref_file=gzip.open(sys.argv[2].replace('chr1','chr'+chr1),'rt')
			print('chr'+chr1)
			ref_line=ref_file.readline().strip()
		ref_indi=ref_line.split('\t')
		refchr=ref_indi[0]
		refpos=int(ref_indi[1])
		if chr1 != refchr:
			print('chr mismatch')
			sys.exit()
		if pos1 > refpos:
			ref_line=ref_file.readline().strip()
		elif pos1 == refpos:
			pinfo=ref_indi[altidx+3]
			bldp=int(ref_indi[7].split(';')[0])
			blrc=int(ref_indi[7].split(';')[altidx+1])
			if bldp ==0:
				blvaf='NA'
			else:
				blvaf=str(round(blrc*100/float(bldp),2))+'%'
			out_file.write(in_line+'\t'+pinfo+'\t'+str(bldp)+';'+str(blrc)+';'+blvaf+'\n')
			in_line=in_file.readline().strip()
		elif pos1 < refpos:
			out_file.write(in_line+'\tNA\tNA\n')
			in_line=in_file.readline().strip()
		prv_chr=chr1
