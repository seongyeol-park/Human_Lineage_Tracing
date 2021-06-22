#This script was made 2018-06-01
#2018-10-08 correction for gzip input
#2021-03-13 make it compatible with python3
import sys,gzip
print(sys.argv[1])
if sys.argv[1][-2:]=='gz':
	in_file=gzip.open(sys.argv[1],'rt')
else:
	in_file=open(sys.argv[1])
sampleid=sys.argv[2]
snp_file=open(sampleid+'.snv.vcf','w')
indel_file=open(sampleid+'.indel.vcf','w')
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		snp_file.write(in_line+'\n')
		indel_file.write(in_line+'\n')
	else:
		in_indi=in_line.split('\t')
		refnt=in_indi[3].split(',')[0]
		altnt=in_indi[4].split(',')[0]
		if len(refnt)== len(altnt):
			snp_file.write(in_line+'\n')
		else:
			indel_file.write(in_line+'\n')
	in_line=in_file.readline().strip()

