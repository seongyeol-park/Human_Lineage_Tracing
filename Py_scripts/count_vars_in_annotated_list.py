import sys
ncol_loca = 11
ncol_exon = 14

ncol_loca = ncol_loca -1
ncol_exon = ncol_exon -1

fn_file=open(sys.argv[1]) # id \t path
out_file = open('variant_count_result.txt','w')
fn_line = fn_file.readline().strip()
out_file.write('sample_id\tn_pointmt\tn_snv\tn_indel\tn_ins\tn_del\tn_exonic\tn_intronic\tn_intergenic\tn_ncRNA\tn_splicing\tn_synony\tn_nonsynony\tn_fsdel\tn_fsins\tn_nfsdel\tn_nfsins\tn_stopgain\tn_stoploss\tn_ukn\n')
while fn_line:
	print(fn_line)
	fn_indi = fn_line.split('\t')
	sample_id = fn_indi[0]
	file_path = fn_indi[1]
	in_file=open(file_path)
	in_line = in_file.readline().strip()
	n_pointmt = 0;n_snv = 0; n_indel = 0; n_ins = 0 ; n_del = 0
	n_exonic = 0 ; n_intronic =0; n_intergenic = 0; n_ncRNA = 0; n_splicing = 0; n_UTR3 = 0; n_UTR5 =0
	n_synony =0; n_nonsynony = 0; n_fsdel =0; n_fsins = 0; n_nfsdel=0; n_nfsins=0; n_stopgain = 0; n_stoploss = 0; n_ukn = 0 
	while in_line:
		if in_line[0]=='#':
			'blank'
		else:
			n_pointmt +=1
			in_indi=in_line.split('\t')
			refnt = in_indi[3]
			altnt = in_indi[4]
			v_loca = in_indi[ncol_loca]
			exon_func = in_indincol_exon]

			if v_loca[0:4]=='exon':
				n_exonic +=1
				if exon_func[0:4] == 'syno':
					n_synony +=1
				elif exon_func[0:4] == 'nons':
					n_nonsynony +=1
				elif exon_func == 'frameshift deletion':
					n_fsdel +=1
				elif exon_func == 'frameshift insertion':
					n_fsins +=1
				elif exon_func == 'nonframeshift deletion':
					n_nfsdel +=1
				elif exon_func == 'nonframeshift insertion':
					n_nfsins +=1
				elif exon_func == 'stopgain':
					n_stopgain +=1
				elif exon_func == 'stoploss':
					n_stoploss +=1
				elif exon_func == 'unknown':
					n_ukn +=1
				else:
					print(f'Exceptional exonic function : {exon_func}')
			elif v_loca[0:4]=='intr':
				n_intronic +=1
			elif v_loca[0:4] == 'inte':
				n_intergenic +=1
			elif v_loca[0:4]=='down':
				n_intergenic +=1
			elif v_loca[0:4] == 'upst':
				n_intergenic +=1
			elif v_loca[0:4] == 'ncRN':
				n_ncRNA +=1
			elif v_loca[0:4] == 'UTR3':
				n_UTR3 +=1
			elif v_loca[0:4] == 'UTR5':
				n_UTR5 +=1
			elif v_loca[0:4] == 'spli':
				n_splicing +=1
			else:
				print(f'Exceptional location : {v_loca}')

			if len(refnt) == 1 and len(altnt) == 1:
				n_snv +=1
			elif len(refnt) > len(altnt) :
				n_del +=1; n_indel +=1
			elif len(refnt) < len(altnt) :
				n_ins +=1; n_indel +=1
			else:
				print('Exeptional variant type.')
				print(in_line)
		in_line = in_file.readline().strip()
	info_list = [sample_id, str(n_pointmt), str(n_snv), str(n_indel), str(n_ins), str(n_del), str(n_exonic), str(n_intronic), str(n_intergenic), str(n_ncRNA), str(n_splicing), str(n_synony), str(n_nonsynony), str(n_fsdel), str(n_fsins), str(n_nfsdel), str(n_nfsins), str(n_stopgain), str(n_stoploss), str(n_ukn)]
	out_file.write('\t'.join(info_list)+'\n')
	fn_line = fn_file.readline().strip()
