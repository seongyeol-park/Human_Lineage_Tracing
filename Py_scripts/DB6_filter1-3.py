#181106: rasm >=2 or pnn >=3 filter out
#181127: column_number_editing; location =='NA' passing as false; only include chromosomes in chr_list; varread < 3 filter out
import sys
print(sys.argv[1])
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.fi3','w')
readlength=int(sys.argv[2])

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
		ltloca=in_indi[14]
		rtloca=in_indi[15]
		clipinfo=in_indi[16]
		insinfo=in_indi[17]
		delinfo=in_indi[18]
		varmm=in_indi[20]
		pninfo=in_indi[21]
		rainfo=in_indi[22]
		p1info=in_indi[23]
		p2info=in_indi[24]
		okg_all=in_indi[30]
		okg_eas=in_indi[31]
		exac_all=in_indi[32]
		exac_eas=in_indi[35]
	
		refread=int(refread)
		varread=int(varread)

		if varread <3 or chr1 not in chr_list:    #filter out
			'blank'
		else:
			locaco300=round(float(readlength)/300**(float(1)/varread),0)  #probability set if 200, then random probability  < 1/200 as error
			if refmq == 'NA':
				refmq=60
			else:
				refmq=float(refmq)
			varmq=float(varmq)
			varmm=float(varmm)
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
				p1vaf=float(p1info.split(';')[6])
			if p2info == 'NA':
				p2vaf=0
			else:
				p2vaf=float(p2info.split(';')[6])
			okg_all=float(okg_all.replace('.','0'))
			okg_eas=float(okg_eas.replace('.','0'))
			exac_all=float(exac_all.replace('.','0'))
			exac_eas=float(exac_eas.replace('.','0'))

                        if refmq-varmq <10 and varmq >= 40 and refmq >=40 and varmm <5 and max(vclipv, vinsv, vdelv) < 90 and p1vaf < 1 and p2vaf < 1: 
				m +=1
				out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()
print(n)
print(m)
