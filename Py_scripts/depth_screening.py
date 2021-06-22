#arg1 bam

target_pos='/home/users/sypark/00_Project/04_QC/01_hg19/04_MPUTargetSelection/02_filter/03_cds_snp_0.49_0.51/hg19_aeg.sites.2015_08_chrall_snp_cds_cut_0.49_0.51.txt.segdup.rmadj.rmBias'

import sys, pysam, numpy

bam_file=pysam.AlignmentFile(sys.argv[1],'rb')
t_file=open(target_pos)
out_file=open(sys.argv[1].replace('.bam','.depthscreen'),'w')
header_list=['#CHROM','POS','depth']
out_file.write('\t'.join(header_list)+'\n')
dp_list=[]
t_line=t_file.readline().strip()
while t_line:
	t_indi=t_line.split('\t')
	chr1=t_indi[0]
	pos1=int(t_indi[1])
	ctarray=bam_file.count_coverage(chr1, pos1-1, pos1, quality_threshold=20) #default: skip unmap, secondary, qcfail, dup
	dp=int(ctarray[0][0])+ctarray[1][0]+ctarray[2][0]+ctarray[3][0]
	dp_list.append(dp)
	out_list=[chr1,str(pos1),str(dp)]
	out_file.write('\t'.join(out_list)+'\n')
	t_line=t_file.readline().strip()

sum_list=[sys.argv[1],str(round(numpy.median(dp_list),2)),str(round(numpy.mean(dp_list),2))]
out_file.write('##input,median,,mean:\t'+'\t'.join(sum_list)+'\n')
