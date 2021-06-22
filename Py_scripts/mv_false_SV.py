import sys, os

os.system('mkdir -p excluded')

f_file=open('/home/users/sypark/00_Project/06_LineageTracing/meta_data/false_SV_list.txt')
f_line=f_file.readline().strip()
false_list=[]
while f_line:
	false_list.append(f_line)
	f_line = f_file.readline().strip()



fn_file=os.popen('ls |grep '+sys.argv[1])
fn_line=fn_file.readline().strip()
while fn_line:
	sv_id = '_'.join((fn_line.split('.')[1]).split('_')[0:5])
	if sv_id in false_list:
		print(fn_line)
		os.system('mv '+fn_line+' excluded/')
	fn_line=fn_file.readline().strip()

