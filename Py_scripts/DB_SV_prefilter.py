#191009 copy from DB3
#200306 modify in short DEL: t_dv < 5 and t_rv == 0 -> t_rv ==0
#200412 NotPASS & n_DV + nRV >= 2 -> filter out / (t_DV + nRV) - ( n_DV + n_RV) < 3 -> filter out

import sys
print(sys.argv[1])
shdel_co=1000


in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.pre_fi','w')
in_line=in_file.readline().strip()
n=0;m=0
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:
		n+=1
		in_indi=in_line.split('\t')
		np_info=in_indi[11]
		np_var=int(np_info.split(';')[0])
		np_n=int(np_info.split(';')[1])
		svtype=in_indi[2][0:3]
		pos1=int(in_indi[1])
		pos2=int((in_indi[7].split('END=')[1]).split(';')[0])
		filter_info = in_indi[6]
		t_dv=int(in_indi[9].split(':')[-3])
		t_rv=int(in_indi[9].split(':')[-1])
		n_dv=int(in_indi[10].split(':')[-3])
		n_rv=int(in_indi[10].split(':')[-1])

		if np_n >=2 :
			'blank'
		elif filter_info != 'PASS' and n_dv + n_rv >=2:
			'blank'
		elif (t_dv + t_rv) - (n_dv + n_rv) < 3:
			'blank'
		elif svtype == 'DEL' and pos2-pos1 < shdel_co and t_rv ==0:
			'blank'
		elif t_dv <3 and t_rv == 0:
			'blank'
		else:	
			m+=1
			out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()
print(n)
print(m)

		
