#191017 written
#200415 subprocess.Popen -> subprocess.getoutput

import subprocess
res_dic={}
total_sigs=[]
info_list=subprocess.getoutput('ls |grep INFO.tsv$')
for line in info_list.split('\n'):
	sampleid=line.split('.INFO.tsv')[0]
	if sampleid not in res_dic.keys():
		res_dic[sampleid]={}
	info_file=open(line)
	info_line=info_file.readline().strip()
	cs=round(float(info_line.split('similarity: ')[1]),2)
	res_dic[sampleid]['cs']=cs

sig_list=subprocess.getoutput('ls |grep signature_exposures.tsv$')
for line in sig_list.split('\n'):
	sampleid = line.split('.signature_exposures')[0]
	res_dic[sampleid]['proportion']={}
	sig_file=open(line)
	sig_line=sig_file.readline().strip()
	sig_line=sig_file.readline().strip() #pass 1st row
	while sig_line:
		sig_indi=sig_line.split('\t')
		signature=sig_indi[0].replace(' ','')
		proportion=float(sig_indi[2])
		total_sigs.append(signature)
		res_dic[sampleid]['proportion'][signature]=proportion
		sig_line=sig_file.readline().strip()

spec_list=subprocess.getoutput('ls |grep signature_spectra.tsv$')
for line in spec_list.split('\n'):
	sampleid = line.split('.signature_spectra')[0]
	spec_file=open(line)
	spec_line=spec_file.readline().strip()
	spec_line=spec_file.readline().strip() # pass 1st row
	n_subs=0
	while spec_line:
		spec_indi=spec_line.split('\t')
		sub_ct=int(spec_indi[2])
		n_subs += sub_ct
		spec_line=spec_file.readline().strip()
	res_dic[sampleid]['n_subs']=n_subs


total_sigs = list(set(total_sigs))
out_file=open('merged_sig_result.txt','w')
header_list=['#sample_id','n_subs','CS']+total_sigs
out_file.write('\t'.join(header_list)+'\n')
for sample in res_dic.keys():
	out_file.write(sample+'\t'+str(res_dic[sample]['n_subs'])+'\t'+str(res_dic[sample]['cs']))
	for sig in total_sigs:
		if sig not in res_dic[sample]['proportion'].keys():
			out_file.write('\t0')
		else:
			out_file.write('\t'+str(res_dic[sample]['proportion'][sig]))
	out_file.write('\n')

	
	
