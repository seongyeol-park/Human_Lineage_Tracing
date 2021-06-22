#!/home/users/sypark/anaconda3/bin/python

#updates
#200307 edited from read_local_reassembly_v2_200304.py
#200406 When I run with multicores, occasional false output was observed. So I removed the option for multicores

import sys,pysam,gzip,collections
from numpy import median


def run_src(in_line, run_mode):
	lt_seq_size=5 # for indel
	snv_seq_size=20 # for snv
	rt_min=5 # Rt minmum sequence for indel
	r_ser=160 # range of searching reference
	max_repeat_unit_seq = 15

	new_lt_size=lt_seq_size
	in_indi=in_line.split('\t')
	chr1=in_indi[0];pos1=int(in_indi[1])
	ref_nt=in_indi[3]
	alt_nt=(in_indi[4].split(',')[0]).split('/')[0]
	called_id = in_indi[args.id_col -1].split(';')[0]
	seq_dic={}
	if len(ref_nt) == len(alt_nt) and len(ref_nt)==1:
		mttype='snv'
		ref_nt=ref_nt[0]
		alt_nt=alt_nt[0]
		new_lt_size=snv_seq_size
		new_rt_size=snv_seq_size
	elif len(ref_nt) == len(alt_nt) and len(ref_nt) > 1:
		mttype=='mnv'
		new_lt_size=snv_seq_size+len(ref_nt)
		new_rt_size=snv_seq_size+len(ref_nt)
	elif len(ref_nt) >  len(alt_nt) and len(alt_nt) == 1:
		mttype='del'
		cor_len=len(alt_nt)-1
		alt_nt=alt_nt[0]
		ref_nt=ref_nt[0:len(ref_nt)-cor_len]
		del_len=len(ref_nt)-1
		del_seq=ref_nt[1:]
		var_seq=ref_nt[1:]
	elif len(ref_nt) < len(alt_nt) and len(ref_nt) == 1:
		mttype='ins'
		cor_len=len(ref_nt)-1
		ref_nt=ref_nt[0]
		alt_nt=alt_nt[0:len(alt_nt)-cor_len]
		ins_len=len(alt_nt)-1
		ins_seq=alt_nt[1:]
		var_seq=alt_nt[1:]
	elif len(ref_nt) != len(alt_nt) and min(len(ref_nt), len(alt_nt)) > 1:
		mttype='cindel'
		new_lt_size=max(len(ref_nt), len(alt_nt))+1
		new_rt_size=max(len(ref_nt), len(alt_nt))+1
	else:
		print('ERROR: unknown mutation type. exiting')
		print(in_line)
		sys.exit()
	


	if mttype =='del' or mttype=='ins':
		ref_seq=r_file.fetch(chr1, pos1-1,pos1+r_ser) 
		repeat_seq=''
		####repeat sequence decomposition
		if len(var_seq) > 1:
			for k in range(1, min(len(var_seq),max_repeat_unit_seq)):
				if len(var_seq) % k==0:
					test_unit=var_seq[0:k]
					for i in range(0,len(var_seq)//len(test_unit)):
						if var_seq[i*k:(i+1)*k]!=test_unit:
							complete_repeat='no'
							break
						else:
							complete_repeat='yes'
				if complete_repeat=='yes':
					repeat_seq=test_unit
					break
		if repeat_seq=='':
			repeat_seq=var_seq

		min_ref_len=0
		for i in range(1, r_ser):
			if ref_seq[i] == repeat_seq[(i-1)%len(repeat_seq)]:
				min_ref_len +=1
			else:
				break
		if mttype== 'ins':
			new_rt_size=max(min_ref_len+len(var_seq)+3, rt_min)
		elif mttype== 'del':
			new_rt_size=max(min_ref_len-len(var_seq)+3, rt_min)
	for i in range(new_lt_size*(-1), new_rt_size+1):
		seq_dic[i]=[]
	ref_n_list=[];var_n_list=[];ukn_n_list=[]
	t_file=pysam.AlignmentFile(in_path_dic[called_id],'rb') #bam file
	for read in t_file.fetch(chr1,pos1-1,pos1):
		if read.is_unmapped == True or read.is_duplicate == True: continue
		var_read='off'
		est_dist=int(in_indi[1])-read.reference_start-1
		rlength=read.infer_query_length()
		if read.cigartuples==None:
			continue
		cigar_list=read.cigartuples
		current_m=0;current_i=0;current_d=0;target_del_stat=0
		for cigar in cigar_list:   #loop for calculate real distance 
			if cigar[0]==0 and (current_m + current_d)  <=est_dist:
				current_m=current_m+cigar[1]
			elif cigar[0]==1 and (current_m + current_d) <=est_dist:
				current_i=current_i+cigar[1]
			elif cigar[0]==2 and (current_m + current_d) <=est_dist:
				if current_m+current_d+cigar[1] > est_dist:
					target_del_stat=1
					break
				else:
					current_d=current_d+cigar[1]
			elif current_m + current_d > est_dist:
				break
			else:
				'blank'
		if target_del_stat==1:
			ref_n_list.append(read.query_name)
		rel_dist=est_dist+current_i-current_d  # start with 0
		
		cigar_i=0;cigar_d=0;cigar_s=0;cigar_h=0; current_m=0;cigar_md=0
		for cigar in cigar_list:   #loop for check presence of clipping, insertion, deletion
			if cigar[0]==0:
				current_m += cigar[1]
			elif cigar[0]==1:
				if mttype == 'ins' and pos1-read.reference_start==current_m+cigar_md and cigar[1]==ins_len and read.query_alignment_sequence[rel_dist+1:rel_dist+1+ins_len]==ins_seq:
					var_read='on'
				else:
					cigar_i=cigar_i+1
			elif cigar[0]==2:
				if mttype == 'del' and pos1-read.reference_start==current_m+cigar_md and cigar[1] == del_len:
					var_read='on'
				else:
					cigar_d=cigar_d+1
				if current_m > 0:
					cigar_md += cigar[1]
			elif cigar[0]==4:
				cigar_s=cigar_s+1
			elif cigar[0]==5:
				cigar_h=cigar_h+1
		if target_del_stat!=1:
			if mttype in ['snv','mnv','cindel'] and read.query_alignment_sequence[rel_dist:rel_dist+len(alt_nt)]==alt_nt:   #var_read
				var_read='on'

		if var_read == 'on':
			var_n_list.append(read.query_name)
			for n in range(0, len(read.query_alignment_sequence)):
				if (n-rel_dist) < new_lt_size*(-1) or n-rel_dist > new_rt_size: continue
				seq_dic[n-rel_dist].append(read.query_alignment_sequence[n])
	prv_len=0;tmp_lt_size=''
	for i in range(new_lt_size*-1,1):  # dictionary check
		if len(seq_dic[i]) == 0 and prv_len !=0:
			print('Sequence Dictionary Error')
			sys.exit()
		if len(seq_dic[i]) !=0 and prv_len ==0:
			tmp_lt_size=abs(i)
		prv_len=len(seq_dic[i])
	prv_len=1;tmp_rt_size=''
	for i in range(0,new_rt_size+1):  # dictionary check
		if len(seq_dic[i]) != 0 and prv_len ==0:
			print('Sequence Dictionary Error')
			sys.exit()
		if len(seq_dic[i])==0:
			tmp_rt_size=i-1
		prv_len=len(seq_dic[i])
	if tmp_lt_size != '':
		new_lt_size=tmp_lt_size
	if tmp_rt_size != '':
		new_rt_size=tmp_rt_size

	cons_seq=''
	for i in range(new_lt_size*-1, new_rt_size+1):
		if len(seq_dic[i])==0:continue
		lt=collections.Counter(seq_dic[i]).most_common(1)[0][0]
		cons_seq=cons_seq+lt
	ori_seq = cons_seq[0:new_lt_size]+ref_nt+cons_seq[new_lt_size+len(alt_nt):]
	if cons_seq=='':
		print('consensus sequence was not made. line skipping')
		print(in_line)
		return('No_cons')
	tcons_list=[];tref_list=[];tukn_list=[]
	for read in t_file.fetch(chr1,pos1-1,pos1):
		if read.is_unmapped == True or read.is_duplicate == True: continue
		if cons_seq in read.query_sequence:
			tcons_list.append(read.query_name)
		elif ori_seq in read.query_sequence:
			tref_list.append(read.query_name)
		else:
			tukn_list.append(read.query_name)
	t_file.close()
	tcons=len(list(set(tcons_list)))
	tref=len(list(set(tref_list)))
	tukn=len(list(set(tukn_list)))
	if tcons + tref == 0:
		tvaf='NA'
	else:
		tvaf=tcons/float(tcons + tref)

	result_list=[]
	for sid in out_path_dic.keys():
		n_file=pysam.AlignmentFile(out_path_dic[sid],'rb') #bam file
		ncons_list=[];nref_list=[];nukn_list=[]
		for read in n_file.fetch(chr1,pos1-1,pos1):
			if read.is_unmapped == True or read.is_duplicate == True: continue
			if cons_seq in read.query_sequence:
				ncons_list.append(read.query_name)
			elif ori_seq in read.query_sequence:
				nref_list.append(read.query_name)
			else:
				nukn_list.append(read.query_name)
		nref=len(list(set(nref_list)))
		ncons=len(list(set(ncons_list)))
		nukn=len(list(set(nukn_list)))
		if ncons + nref < 5:
			nvaf = 'NA'
		else:
			nvaf = round(ncons*100/float(ncons + nref),2)
		if run_mode == 'vafonly':
			result_list.append(str(nvaf))
		else:
			result_list.append(';'.join([str(nref), str(ncons), str(nukn), str(nvaf)]))
		n_file.close()
	if mttype == 'snv':
		repeat_seq='.'
	info_list=[in_line, called_id]+result_list
	return info_list


#collect line
if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(
		description = "This script counts the number of vairant reads, reference reads, and unknown reads using read local ressembly from multiple BAM files. This script was modified to get multiple input and output BAMs for merged VCFs.")
	parser.add_argument("vcf", type=str, help="input vcf file")
	parser.add_argument("id_col", type=int, help="column number containing called ID")
	parser.add_argument("input_id_bam", type=str, help="inputs which are consistent with id_col (id1 \t bam_path1\n id2 \t bam_path2 \n ...)")
	parser.add_argument("output_id_bam", type=str, help="samples you want information of (id1 \t bam_path1 \n id2 \t bam_path2 \n ...)")
	parser.add_argument("-r","--reference", type=str, default='/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta', help="default[/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta]")
	parser.add_argument("-v","--vafonly", action = "store_true", help = "Write only VAF information rather than read counts")
	args=parser.parse_args()
	
	print(sys.argv[0])
	print(args.vcf)
	if args.vcf[-3:]=='.gz':
		in_file=gzip.open(args.vcf)
	else:
		in_file=open(args.vcf)
	r_file=pysam.FastaFile(args.reference)
	if args.vcf[-3:]=='.gz':
		out_file=open(args.vcf[:-3]+'.mvaf_test','w')
	else:
		out_file=open(args.vcf+'.mvaf_test','w')
	id_file=open(args.input_id_bam)
	id_line=id_file.readline().strip()
	in_path_dic={}
	while id_line:
		id_indi=id_line.split('\t')
		sample_id = id_indi[0]
		sample_path = id_indi[1]
		in_path_dic[sample_id]=sample_path
		id_line = id_file.readline().strip()
	id_file.close()
	
	id_file=open(args.output_id_bam)
	id_line=id_file.readline().strip()
	out_path_dic={}
	while id_line:
		id_indi=id_line.split('\t')
		sample_id = id_indi[0]
		sample_path = id_indi[1]
		out_path_dic[sample_id]=sample_path
		id_line = id_file.readline().strip()
	id_file.close()

	if args.vafonly:
		run_mode = "vafonly"
	else:
		run_mode = "full"


	line_list=[]
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0:4]=='#CHR':
			out_file.write(in_line+'\tcalled_sample')
			if run_mode == 'vafonly':
				for sid in out_path_dic.keys():
					out_file.write('\t'+sid+'_vafpct')
			else:
				for sid in out_path_dic.keys():
					out_file.write('\t'+sid+'_ref;'+sid+'_var;'+sid+'_ukn;'+sid+'_vafpct')
			out_file.write('\n')
		elif in_line[0]=='#':
			out_file.write(in_line+'\n')
		else:   
			res_list = run_src(in_line, run_mode)
			if res_list == 'No_cons':
				'blank'
			else:
				out_file.write('\t'.join(res_list)+'\n')
		in_line=in_file.readline().strip()
	
	print('done')
