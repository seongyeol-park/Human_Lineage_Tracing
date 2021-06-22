#2019-12-20 copied and edited from DB5_mpu_snv_sampleN_Dpmore10.py

#You should update sample_name_list, exp_list, blood_list

import sys,gzip
in_file=gzip.open(sys.argv[1])  #*.snv.gz


# chr \t pos \t ref \t A_readN;sampleN;sampleNR2;sampleNR3;sampleNR3_names \t G... \t C... \t T...\t blood_DP;A;G;C;T n
in_line=in_file.readline().strip()
nt_list=['A','G','C','T']
sample_name_list=['10_ALL13-4_006H8','10_ALL13-8_001G5','10_ALL17-3_001D3','10_ALL17-5_002B3','10_ALL19-2_001B5','10_ARL10-4_001D4-B6','10_ARL10-4_001D4','10_ARL12-1_001B3','10_ARL12-2_001C12','10_ARL12-4_001H11','10_ARL1-3_003G11','10_ARL13-1_001G8','10_ARL1-5_003B3','10_ARL23-6_002E11','10_ARL3-2_006B8','10_ARL6-1_001A9','10_ARL6-2_001E7','10_ARL7-2_005B2','10_ARL8-4_001D12','10_ARL9-2_005G8','10_ARL9-3_001D5','10_ARLFb2-5_002G8','10_Blood_3s_merged','10_LA14-2_001E7','10_LA14-2_002D1-A5','10_LA14-2_002D1','10_LA15-2_001D7','10_LAFb11-3_001G8','10_LAFb16-2_001D8','10_LAFb16-3_001E11','10_LAFb8-2_001D4','10_PLL2-2_001F4','10_PLL4-4_003F8','10_PLL6-3_002A11','10_PLL6-3_002D6','10_PRL10-2_001G10','10_PRL10-3_001G7','10_PRL10-5_001E2','10_PRL12-5_001C5','10_PRL16-1_001E6','10_PRL2-4_002C4','10_PRL3-1_001B8-E10','10_PRL3-1_001B8','10_PRL3-1_001E4','10_PRL3-2_002E11','10_PRL3-3_001B11','10_PRL4-1_004C1','10_PRL4-4_001E8','10_PRLM2-3_001D10','10_RA9-3_005G11','10_RA9-4_005B8','10_UO3-2_001C6']
total_n=len(sample_name_list)
sample_list=range(1,total_n+1)  # index of clonally expanded samples
sample_list=[x+2 for x in sample_list]
exp_list=[22]  #index of excluded samples
exp_list=[x+2 for x in exp_list]
blood_list=[23]  # index of blood samples
blood_list=[x+2 for x in blood_list]
sample_list = list(set(sample_list) - set(exp_list) - set(blood_list))
filtered_n=len(sample_list)
out_file=gzip.open(sys.argv[1].replace('.gz','.'+str(filtered_n)+'N.gz'),'w')
while in_line:
	in_indi=in_line.split('\t')
	bldp=0;blA=0;blG=0;blC=0;blT=0
	out_file.write(in_indi[0]+'\t'+in_indi[1]+'\t'+in_indi[2])
	for nt in nt_list:
		snr1=0
		snr2=0
		snr3=0
		snr3_name='None'
		sread=0
		for i in sample_list:
			s_name=sample_name_list[i-3]
			if nt+'=' in in_indi[i]:
				readn=int((in_indi[i].split(nt+'=')[1]).split(';')[0])
				snr1+=1
				sread += readn
				if readn >=2:
					snr2 +=1
				if readn >=3:
					snr3 +=1
					if snr3_name=='None':
						snr3_name=s_name
					else:
						snr3_name=snr3_name+';'+s_name
		if snr3==filtered_n:
			snr3_name='All'
		info_list=[str(sread), str(snr1), str(snr2), str(snr3), snr3_name]
		out_file.write('\t'+';'.join(info_list))
	for i in blood_list:
		dp=int((in_indi[i].split('dp=')[1]).split(';')[0])
		bldp += dp
		if 'A=' in in_indi[i]:
			Ar=int((in_indi[i].split('A=')[1]).split(';')[0])
			blA += Ar
		if 'G=' in in_indi[i]:
			Gr=int((in_indi[i].split('G=')[1]).split(';')[0])
			blG += Gr
		if 'C=' in in_indi[i]:
			Cr=int((in_indi[i].split('C=')[1]).split(';')[0])
			blC += Cr
		if 'T=' in in_indi[i]:
			Tr=int((in_indi[i].split('T=')[1]).split(';')[0])
			blT += Tr
	binfo_list=[str(bldp),str(blA),str(blG),str(blC),str(blT)]
	out_file.write('\t'+';'.join(binfo_list)+'\n')
	in_line=in_file.readline().strip()
