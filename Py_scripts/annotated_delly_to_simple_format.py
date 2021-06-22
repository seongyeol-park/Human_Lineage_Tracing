#20190930 written
#20191010 add varn refn

#arg1: input delly vcf

import sys, os
tcol=10
ncol=11

tcol=tcol-1
ncol=ncol-1

fn_line=sys.argv[1]
in_file=open(fn_line)
out_file=open(fn_line+'.simple_edit','w')
header_list=['#CHR1','POS1','CHR2','POS2','mh','terinfo','svtype', 'dellyid','ref1','ref2','total_var','split_var','SA_var','vaf1','vaf2']
out_file.write('\t'.join(header_list)+'\n')

in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		#chr1=in_indi[0]
		#pos1=int(in_indi[1])
		#svtype=in_indi[2][0:3]
		delly_id=in_indi[2]
		filter_info=in_indi[6]
		precise_info=in_indi[7].split(';')[0]
		mapq=int((in_indi[7].split('MAPQ=')[1]).split(';')[0])
		#chr2=(in_indi[7].split('CHR2=')[1]).split(';')[0]
		#pos2=int((in_indi[7].split('END=')[1]).split(';')[0])
		#ori=(in_indi[7].split('CT=')[1]).split(';')[0]
		tDV=int(in_indi[tcol].split(':')[9])
		tRV=int(in_indi[tcol].split(':')[11])
		nDV=int(in_indi[ncol].split(':')[9])
		nRV=int(in_indi[ncol].split(':')[11])
		chr1=in_indi[14]
		pos1=int(in_indi[15])
		chr2=in_indi[16]
		pos2=int(in_indi[17])
		mh=in_indi[18]
		ori=in_indi[19]
		svtype=in_indi[20]			
		frag_info=in_indi[26]
	
		info_list=[chr1, str(pos1), chr2, str(pos2), mh, ori, svtype, delly_id] + frag_info.split(';')
		rvs_info_list=[chr2, str(pos2), chr1, str(pos1), ori[-1]+'to'+ori[0], svtype, delly_id] + frag_info.split(';')
		chr1n=int(chr1.replace('X','23').replace('Y','24'))
		chr2n=int(chr2.replace('X','23').replace('Y','24'))
		if chr1n > chr2n:
			out_file.write('\t'.join(rvs_info_list)+'\n')
		else:
			out_file.write('\t'.join(info_list)+'\n')
	in_line=in_file.readline().strip()
