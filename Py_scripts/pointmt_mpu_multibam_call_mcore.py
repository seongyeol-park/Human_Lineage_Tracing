#args1: input
#args2: bam_path_list
#args3: id of bams
#args5: number of cores for multiprocess

#191026 copied and edited from DB6 filter3 file, bam list as sys.argv[2], run_serial calling and editing, lowDp NA
#200525 Major error correction

import sys,os,gzip,subprocess
from multiprocessing import Pool


def run_mpu_call(chr1, pos1):
	bam_list =sys.argv[2]
	cmd='samtools mpileup -AB -Q 20 -q 20 -d 50000 -b '+bam_list+' -f /home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta -r '+chr1+':'+pos1+'-'+pos1
	print(cmd)
	res=subprocess.Popen(cmd, stdout = subprocess.PIPE, shell = True)
	res.wait()
	output, error = res.communicate()
	errcode = res.returncode
	if errcode != 0:
		print(error)
		sys.exit(1)
	else:
		return (output)

def run_mpu_call_mp(inputs,ncor):
	pool = Pool(processes = ncor)
	res = pool.starmap(run_mpu_call, inputs)
	return(res)

def get_sampleN(input_file):
	in_file = open(input_file)
	in_line = in_file.readline().strip()
	n=0
	while in_line:
		in_line=in_file.readline().strip()
		n = n+1
	in_file.close()
	return(n)

def annotate_SNV(t_chr,t_pos,refnt,altnt,snv_dic,dic_bin):
	refnt = refnt[0] 
	altnt = altnt[0]
	t_line = snv_dic[t_chr][t_pos//dic_bin][t_pos]
	ref_list = t_line.split('\t')
	out_info_list=[]
	for info in ref_list:
		c_dp=int((info.split('dp=')[1]).split(';')[0])
		if refnt+'=' in info:
			c_ref=int((info.split(refnt+'=')[1]).split(';')[0])
		else:
			c_ref=0
		if altnt+'=' in info:
			c_var=int((info.split(altnt+'=')[1]).split(';')[0])
		else:
			c_var=0
		if c_dp==0:
			c_vaf='NA'
		else:
			c_vaf=str(round(c_var*100/float(c_dp),2))
		c_info_list=[str(c_dp), str(c_ref), str(c_var), str(c_vaf)]
		c_info='\t'.join(c_info_list)
		out_info_list.append(c_info)
	return '\t'.join(out_info_list)

def annotate_indel(t_chr,t_pos,refnt,altnt,indel_dic,dic_bin):
	out_info_list=[]
	if len(refnt) > len(altnt):
		mttype='del'
		del_seq = refnt[1:]
		var_id = '-' + str(len(del_seq))+del_seq
	elif len(refnt) < len(altnt):
		mttype='ins'
		ins_seq = altnt[1:]
		var_id = '+' + str(len(ins_seq))+ins_seq
	
	t_line = indel_dic[t_chr][t_pos//dic_bin][t_pos]
	ref_list = t_line.split('\t')
	out_info_list=[]
	for info in ref_list:
		c_dp=int((info.split('dp=')[1]).split(';')[0])
		if var_id in info:
			c_var = int(info.split(var_id+'(')[1].split(')')[0])
		else:
			c_var = 0
		c_ref = c_dp - c_var
		if c_dp==0:
			c_vaf='NA'
		else:
			c_vaf=str(round(c_var*100/float(c_dp),2))
		c_info_list=[str(c_dp), str(c_ref), str(c_var), str(c_vaf)]
		c_info='\t'.join(c_info_list)
		out_info_list.append(c_info)
	return '\t'.join(out_info_list)

if __name__ == '__main__':

	ref_fai='/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fa.fai'
	dic_bin=1000000

	pos_file=open(sys.argv[1])
	bam_list = sys.argv[2]
	id_list = sys.argv[3]
	ncor = int(sys.argv[4])
	sampleN = get_sampleN(bam_list)
	print(f"No. of samples : {sampleN}")
	snv_input_list=[]
	indel_input_list=[]
	print('classify SNV and indel')
	pos_line=pos_file.readline().strip()
	while pos_line:
		if pos_line[0]=='#':
			'blank'
		else:
			pos_indi=pos_line.split('\t')
			chr1=pos_indi[0]
			pos1=pos_indi[1]
			refnt = pos_indi[3]
			altnt = pos_indi[4]
			if len(refnt) ==  len(altnt):
				snv_input_list.append((chr1, pos1))
			else:
				indel_input_list.append((chr1,pos1))
		pos_line=pos_file.readline().strip()

	print('run SNV MPU')
	out_file=open(sys.argv[1]+'.snv.mpu', 'w')
	res_lines = run_mpu_call_mp(snv_input_list, ncor)
	for res_line in res_lines:
		txt_out = res_line.decode('utf-8').rstrip()
		out_file.write(txt_out+'\n')
	out_file.close()

	print('run indel MPU')
	out_file=open(sys.argv[1]+'.indel.mpu','w')
	res_lines = run_mpu_call_mp(indel_input_list, ncor)
	for res_line in res_lines:
		txt_out = res_line.decode('utf-8').rstrip()
		out_file.write(txt_out+'\n')
	out_file.close()


	print("Run MPU call")
	cmd1='python /home/users/sypark/01_Python_files/normal_matrix/pileup_calling_snv.py '+sys.argv[1]+'.snv.mpu '+str(sampleN)
	os.system(cmd1)
	cmd2='python /home/users/sypark/01_Python_files/normal_matrix/pileup_calling_indel_190705.py '+sys.argv[1]+'.indel.mpu '+str(sampleN)
	os.system(cmd2)

	print('making SNV and indel dic')
	size_file=open(ref_fai)
	snv_dic={}
	indel_dic={}
	size_line=size_file.readline().strip()
	while size_line:
		size_indi=size_line.split('\t')
		snv_dic[size_indi[0]]={}
		indel_dic[size_indi[0]]={}
		for i in range(0, int(size_indi[1])//dic_bin+1):
			snv_dic[size_indi[0]][i]={}
			indel_dic[size_indi[0]][i]={}
		size_line=size_file.readline().strip()
	
	ref_file=open(sys.argv[1]+'.snv.mpu.snv')
	ref_line=ref_file.readline().strip()
	while ref_line:
		ref_indi=ref_line.split('\t')
		ref_chrom = ref_indi[0]
		ref_pos = int(ref_indi[1])
		snv_dic[ref_chrom][ref_pos//dic_bin][ref_pos]='\t'.join(ref_indi[3:])
		ref_line = ref_file.readline().strip()
	ref_file.close()

	ref_file=open(sys.argv[1]+'.indel.mpu.indel')
	ref_line=ref_file.readline().strip()
	while ref_line:
		ref_indi=ref_line.split('\t')
		ref_chrom = ref_indi[0]
		ref_pos = int(ref_indi[1])
		indel_dic[ref_chrom][ref_pos//dic_bin][ref_pos]='\t'.join(ref_indi[3:])
		ref_line = ref_file.readline().strip()
	ref_file.close()

	
	print('annotating')
	in_file=open(sys.argv[1])
	out_file=open(sys.argv[1]+'.'+str(sampleN)+'sCall','w')
	in_line=in_file.readline().strip()
	n=0
	nt_list=['A','G','C','T']
	id_file=open(id_list)
	id_line=id_file.readline().strip()
	header_list=[]
	while id_line:
		header_list=header_list+[id_line+'_dp',id_line+'_ref',id_line+'_var',id_line+'_vafpct']
		id_line=id_file.readline().strip()
	
	if len(header_list) != sampleN*4:
		print('ERROR: inconsisten sampleN and id_list')
		sys.exit(1)
	
	while in_line:
		in_indi=in_line.split('\t')
		if in_line[0:6]=='#CHROM':
			out_file.write(in_line+'\t'+'\t'.join(header_list)+'\n')
		elif in_line[0]=='#':
			out_file.write(in_line+'\n')
		else:
			t_chr = in_indi[0]
			t_pos = int(in_indi[1])
			refnt = in_indi[3]
			altnt = in_indi[4]
			if len(refnt) == len(altnt) and len(refnt) >1:
				print('WARNING: MNV will treated as SNV with first nucleotide')
			if len(refnt) != len(altnt) and len(refnt) > 1 and len(altnt) >1:
				print('WARNING: No information will be written at complex indel.')
				out_file.write(in_line+'\t'+'\t'.join(['NA']*4*sampleN)+'\n')
				in_line = in_file.readline().strip()
				continue
			if len(refnt) == len(altnt):
				t_info = annotate_SNV(t_chr, t_pos, refnt, altnt,snv_dic, dic_bin)
				out_file.write(in_line+'\t'+t_info+'\n')
			else:
				t_info = annotate_indel(t_chr, t_pos, refnt, altnt,indel_dic, dic_bin)
				out_file.write(in_line+'\t'+t_info+'\n')
		in_line = in_file.readline().strip()
