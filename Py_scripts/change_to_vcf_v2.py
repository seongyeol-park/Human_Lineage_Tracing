#200124 printing error correction: AD <-> RD

import sys
print(sys.argv[1])
in_file=open(sys.argv[1])
sampleid=sys.argv[2]
out_file=open(sampleid+'.edited.vcf','w')
out_file.write('##fileformat=VCFv4.1\n')
out_file.write('##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">\n')
out_file.write('##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">\n')
out_file.write('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">\n')
out_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sampleid+'\n')

in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=in_indi[1]
		ID='.'
		refnt=in_indi[3]
		altnt=in_indi[4]
		qual='.'
		fi='PASS'
		trc=in_indi[5].split(';')[0]
		tac=in_indi[5].split(';')[1]
                info_list=['SOMATIC']
                info=';'.join(info_list)
                format_list=['RD','AD']
                format_info=':'.join(format_list)
                sampleinfo_list=[trc, tac]
                sampleinfo=':'.join(sampleinfo_list)
		final_list=[chr1,pos1,ID,refnt,altnt,qual,fi,info,format_info, sampleinfo]
		out_file.write('\t'.join(final_list)+'\n')
	in_line=in_file.readline().strip()

