#Arg1:input
#Arg2: tumor cell fraction

#Original script: SV_mutCN_CCF_annot.py
#2018-06-14 modification format for SNV
#2018-06-28 chrY error resolve
#2018-08-09 add exception for 'NA' d/t pcawg call
#2019-04-11 resolve logarithm error which occur at purity 1.0
#2020-01-25 edit column number for lineage tracing project
import math
def nCr(n,r):
	f=math.factorial
	return math.log(f(n),10)-math.log(f(r),10)-math.log((f(n-r)),10)
#/  pow(2*math.pi*r,0.5)*pow(r/math.exp(1),r) / pow(2*math.pi*n-r,0.5)*pow(n-r/math.exp(1),n-r)

def calc_scF(tcf, wt_count, var_count, CN, mCN):
	tcf=float(tcf); wt_count=int(wt_count); var_count=int(var_count); CN=float(CN); mCN=float(mCN)
	read_depth=wt_count + var_count
	Rfraction_t=CN*tcf/(CN*tcf+(1-tcf)*2)
	vaf=var_count/float(wt_count+var_count)
	global mutCN, max_mutCN, scF
	if Rfraction_t<=0:
		mutCN=0;max_mutCN=0;scF='.'
	else:
		mutCN=vaf*CN/Rfraction_t

		if CN==2 and mCN==1:
			max_mutCN=1
			scF=mutCN
		else:
			max_binomP=0
			max_mutCN=1
			for this_CN in range(1,max(int(round(CN-mCN,0))+1,1)):
				this_binomP=0
				for readnum in range(var_count,read_depth):
					if Rfraction_t == 1:
						Rfraction_t = 0.999
					this_P1=nCr(read_depth,readnum)+readnum*math.log(Rfraction_t,10)+(read_depth-readnum)*math.log(1-Rfraction_t,10)
					this_P2=nCr(readnum,var_count)+math.log(this_CN/float(CN),10)*var_count+math.log(max(1-this_CN/float(CN),0.0001),10)*(readnum-var_count)
					this_binomP+=pow(10,this_P1+this_P2)
	
				if this_binomP>= max_binomP:
					max_binomP=this_binomP
					max_mutCN=this_CN
			scF=round(mutCN/float(max_mutCN),3)


import sys
fn=sys.argv[1]
print(fn)
inputfile=open(fn)
line=inputfile.readline()

ofn=fn+'.scF' #subclonal Fraction
outputfile=open(ofn,"w")
tcf=float(sys.argv[2])

while line:
	if line[0:6]=='#CHROM':
		outputfile.write(line.rstrip()+"\tdt_mutCN\tdt_mutCH\tdt_CCF\n")
		line=inputfile.readline()
		continue
	
	else:
		line_split=line.rstrip().split("\t")
		wt_count1=int(line_split[6].split(';')[0])
		var_count1=int(line_split[6].split(';')[1])
		if line_split[7]=='.' or line_split[7]=='NA':
			info_list=['.','.','.']
		else:
			seqztCN1=int(line_split[7])
			seqzmCN1=int(line_split[9])
			calc_scF(tcf, wt_count1, var_count1, seqztCN1, seqzmCN1)
			mutCN1=mutCN; max_mutCN1=max_mutCN; scF1=scF
			info_list=[str(mutCN1),str(max_mutCN1), str(scF1)]
		outputfile.write(line.rstrip()+'\t'+'\t'.join(info_list)+'\n')
	line=inputfile.readline()
print ("done")

