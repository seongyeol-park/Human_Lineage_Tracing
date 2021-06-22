import sys, subprocess, argparse,os

py_file='/home/users/sypark/01_Python_files/normal_matrix/AddNpanelToVCFindel_vaf_fast_multi.py'
file_list=['/home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/01_SNU30/snuN30.30s.q0Q0.chr1.mpileup.indel.edit.gz','/home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/04_BGI24/bgiN24.24s.q0Q0.chr1.mpileup.indel.edit.gz','/home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/02_PCAWG36/pcawgN36.36s.q0Q0.chr1.mpileup.indel.edit.gz','/home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/03_SNU24/snuN24.24s.q0Q0.chr1.mpileup.indel.edit.gz','/home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/05_DFCI24/DFCI24.24s.q0Q0.chr1.mpileup.indel.edit.gz']
name_list=['snuN30','bgiN24','pcawgN36','snuN24','dfciN24']


parser = argparse.ArgumentParser(
	description = 'this script annotate indel normal panel to indel vcf')
parser.add_argument("list", help="list of input vcf")
parser.add_argument("-i","--items", default='snuN30;bgiN24;pcawgN36', help="write name of normal panel which you want to annotate. default = 'snuN30;bgiN24;pcawgN36'. possible_list='"+';'.join(name_list)+'\n')
args=parser.parse_args()


current_list=[]
for np in args.items.split(';'):
	print(np)
	out_file=open('tmp_input','w')
	in_file=open(args.list)
	in_line=in_file.readline().strip()
	while in_line:
		if len(current_list) >= 1:
			out_file.write(in_line+'.'+'.'.join(current_list)+'\n')
		else:
			out_file.write(in_line+'\n')
		in_line = in_file.readline().strip()
	in_file.close()
	out_file.close()
	index_num=name_list.index(np)
	cmd_line='python '+py_file + ' tmp_input '+file_list[index_num]+' '+name_list[index_num]
	print(cmd_line)
	res=subprocess.Popen(cmd_line, shell=True)
	if res.wait() != 0:
		output, error = res.communicate()
		print(error)
	current_list.append(np)

os.system('rm -f tmp_input')
