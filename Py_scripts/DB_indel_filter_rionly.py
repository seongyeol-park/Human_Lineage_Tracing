#arg1: input
#arg2: readlength
#arg3: OL column numbers (e.g. "46,47,48")
#arg4: rasm_version (1 or 2)

#181106: rasm >=2 or pnn >=3 filter out
#181127: column_number_editing; location =='NA' passing as false; only include chromosomes in chr_list; varread < 3 filter out
#190817 ncons < 3 only
#190927: other lineage corrected tumor variant reads >= 1;  filter out
#200229: other lineage n_consen >=1 in any two OL samples -> filter out
#200303: other lineage n_conse >=2 -> filter out
#200305: other lineage n_conse >=3 -> filter out
#200306: add max(ol1_nref, ol2_nref) <=1 -> filter out
#200324: get multiple OLcolumns by sys.argv;if (max(info_list) > 1 and refMinMQ <40) or min(info_list) > 1; npanel 1 -> 0.5
#200326: add option of rasm version


import sys
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.ri_fi','w')
readlength=int(sys.argv[2])
OLcols=sys.argv[3].split(',')
OLcols = [int(x)-1 for x in OLcols]
rasm_v = sys.argv[4]

chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
n=m=0
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		out_file.write(in_line+'\n')
	else:	
		n +=1
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		refread=in_indi[10]
		varread=in_indi[11]
		refmq=in_indi[12]
		varmq=in_indi[13]
		ltloca=in_indi[16]
		rtloca=in_indi[17]
		clipinfo=in_indi[18]
		insinfo=in_indi[19]
		delinfo=in_indi[20]
		varmm=in_indi[22]
		pninfo=in_indi[23]
		rainfo=in_indi[24]
		p1info=in_indi[25]
		p2info=in_indi[26]
		okg_all=in_indi[32]
		okg_eas=in_indi[33]
		exac_all=in_indi[34]
		exac_eas=in_indi[37]
	
		refread=int(refread)
		varread=int(varread)

		if varread <3 or chr1 not in chr_list:    #filter out
			'blank'
		else:
			locaco300=round(float(readlength)/300**(float(1)/varread),0)  #probability set if 200, then random probability  < 1/200 as error
			if refmq == 'NA':
				refMedMQ = 60
				refMinMQ = 60
			else:
				refMedMQ = float(refmq.split(';')[1])
				refMinMQ = float(refmq.split(';')[0])
			if varmq == 'NA':
				varMedMQ = 60
			else:
				varMedMQ = float(varmq.split(';')[1])
			varMedMM=float(varmm.split(';')[1])
			if ltloca == 'NA' or rtloca == 'NA':
				in_line=in_file.readline().strip()
				continue
			lmax=int(ltloca.split(';')[2])
			rmax=int(rtloca.split(';')[2])
			rclipv=clipinfo.split(';')[1]
			if rclipv == 'NA':
				rclipv = 0
			else:
				rclipv = float(rclipv)
			vclipv=float(clipinfo.split(';')[3])
			rinsv=insinfo.split(';')[1]
			if rinsv == 'NA':
				rinsv = 0
			else:
				rinsv = float(rinsv)
			vinsv=float(insinfo.split(';')[3])
			rdelv=delinfo.split(';')[1]
			if rdelv == 'NA':
				rdelv = 0
			else:
				rdelv = float(rdelv)
			vdelv=float(delinfo.split(';')[3])
			pnvaf=pninfo.split(';')[2]
			if pnvaf == 'NA':
				pnvaf =0
			else:
				pnvaf = float(pnvaf)
			pnn=int(pninfo.split(';')[1])
			ncons=int(rainfo.split(';')[8])
			if p1info == 'NA':
				p1vaf=0
			else:
				p1vaf=float(p1info.split(';')[4])
			if p2info == 'NA':
				p2vaf=0
			else:
				p2vaf=float(p2info.split(';')[4])
			okg_all=float(okg_all.replace('.','0'))
			okg_eas=float(okg_eas.replace('.','0'))
			exac_all=float(exac_all.replace('.','0'))
			exac_eas=float(exac_eas.replace('.','0'))
			info_list=[]
			if refMedMQ - varMedMQ <10 and varMedMQ >= 40 and refMedMQ >=40 and varMedMM <5 and max(vclipv, vinsv, vdelv) < 90  and ncons < 3:
				m +=1
				out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()
print(n)
print(m)
sys.stdout.flush()
