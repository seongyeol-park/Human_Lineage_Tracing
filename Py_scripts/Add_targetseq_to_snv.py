#Arg1 SNV
#Arg2 mpu snv call
#Arg3 sample  name for header (this should be same order with Arg2)

import sys,gzip
in_file=open(sys.argv[1])  # vcf file
in_line=in_file.readline().strip()
n=0
nt_list=['A','G','C','T']


#count column number
ref_file=open(sys.argv[2])
ref_line=ref_file.readline().strip()
col_num = len(ref_line.split('\t'))-3

#out_file naming
out_file=open(sys.argv[1]+'.'+str(col_num)+'sCall','w')

#make sample name list
sn_file=open(sys.argv[3])
sn_line=sn_file.readline().strip()
sn_list=[]
while sn_line:
	sn_list.append(sn_line)
	sn_line = sn_file.readline().strip()


#annotate
ref_indi=ref_line.split('\t')
while in_line:
	in_indi=in_line.split('\t')
	if in_line[0:6]=='#CHROM':
		out_file.write(in_line+'\t'+'\t'.join(sn_list)+'\n')
		in_line=in_file.readline().strip()
	elif in_line[0]=='#':
		out_file.write(in_line+'\n')
		in_line=in_file.readline().strip()
	elif ref_line=='':
		out_file.write(in_line+'\t'+'\t'.join(['NA']*col_num)+'\n')
		in_line=in_file.readline().strip()
	elif in_indi[3] not in nt_list:
		out_file.write(in_line+'\t'+'\t'.join(['NA']*col_num)+'\n')
		in_line=in_file.readline().strip()
	else:
		ref_indi=ref_line.split('\t')
		chr1=int(in_indi[0].replace('X','23').replace('Y','24'))
		chr2=int(ref_indi[0].replace('X','23').replace('Y','24'))
		refnt=in_indi[3]
		varnt=in_indi[4]
		ref_list=ref_indi[3:]
		if chr1==chr2:
			if int(in_indi[1])==int(ref_indi[1]):
				out_file.write(in_line)
				for info in ref_list:
					c_dp=int((info.split('dp=')[1]).split(';')[0])
					if refnt+'=' in info:
						c_ref=int((info.split(refnt+'=')[1]).split(';')[0])
					else:
						c_ref=0
					if varnt+'=' in info:
						c_var=int((info.split(varnt+'=')[1]).split(';')[0])
					else:
						c_var=0
					if c_dp==0:
						c_vaf='NA'
					else:
						c_vaf=str(round(c_var*100/float(c_dp),2))+'%'
					c_info_list=['dp='+str(c_dp),'ref='+str(c_ref),'var='+str(c_var),'vaf='+c_vaf]
					c_info=';'.join(c_info_list)
					out_file.write('\t'+c_info)
				out_file.write('\n')
				in_line=in_file.readline().strip()
				ref_line=ref_file.readline().strip()
				n+=1
			elif int(in_indi[1])>int(ref_indi[1]):
				ref_line=ref_file.readline().strip()
				n+=1
			elif int(in_indi[1])<int(ref_indi[1]):
				out_file.write(in_line+'\t'+'\t'.join(['NA']*col_num)+'\n')
				in_line=in_file.readline().strip()
		elif chr1>chr2:
			ref_line=ref_file.readline().strip()
			n+=1
		elif chr1<chr2:
			out_file.write(in_line+'\t'+'\t'.join(['NA']*col_num)+'\n')
			in_line=in_file.readline().strip()
	if n!=0 and n%10000000==0:
		print(n)
