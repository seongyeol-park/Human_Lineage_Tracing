#191015 edited: remove precise info
#200420 mkdir -p igv_gff and save output files into the igv_gff folder

import sys,os
in_file=open(sys.argv[1])
os.system('mkdir -p igv_gff')

out_file=open('igv_gff/'+sys.argv[1]+'.gff','w')
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=in_indi[1]
		chr2=in_indi[2]
		pos2=in_indi[3]
		term1=in_indi[5].split('to')[0]
		term2=in_indi[5].split('to')[1]
		svtype=in_indi[6]
		if term1=='3':
			pos1=int(pos1)
			out_file.write(chr1+'\tmarker\t'+svtype+'\t'+str(pos1-99)+'\t'+str(pos1)+'\t.\t+\t.\t'+'ID='+'_'.join(in_indi[0:4])+'_'+in_indi[5]+'_'+svtype+'_1\n')
		elif term1=='5':
			pos1=int(pos1)
			out_file.write(chr1+'\tmarker\t'+svtype+'\t'+str(pos1)+'\t'+str(pos1+99)+'\t.\t-\t.\t'+'ID='+'_'.join(in_indi[0:4])+'_'+in_indi[5]+'_'+svtype+'_1\n')
		if chr2 != '.':			
			if term2=='3':
				pos2=int(pos2)
				out_file.write(chr2+'\tmarker\t'+svtype+'\t'+str(pos2-99)+'\t'+str(pos2)+'\t.\t+\t.\t'+'ID='+'_'.join(in_indi[0:4])+'_'+in_indi[5]+'_'+svtype+'_2\n')
			elif term2=='5':
				pos2=int(pos2)
				out_file.write(chr2+'\tmarker\t'+svtype+'\t'+str(pos2)+'\t'+str(pos2+99)+'\t.\t-\t.\t'+'ID='+'_'.join(in_indi[0:4])+'_'+in_indi[5]+'_'+svtype+'_2\n')
	in_line=in_file.readline().strip()

