#181030 made by SPark
#181127 column error correction when using -s option.
#190315 make this run with hg38
#190318 add refGene of mm10
import sys,os,argparse

possible_list=['refGene','cytoBand','1000g2015aug_all','1000g2015aug_afr','1000g2015aug_amr','1000g2015aug_eas','1000g2015aug_eur','1000g2015aug_sas','exac03','exac03nontcga','esp6500siv2_all','esp6500siv2_aa','esp6500siv2_ea','snp142','avsnp150','clinvar_20180603','dbnsfp35a','dbscsnv11','intervar_20180118','cosmic86_coding','cosmic86_noncoding']
parser = argparse.ArgumentParser()
parser.add_argument("vcf", help="vcf")
parser.add_argument("-i","--items", default='refGene', help="Write items which you want to annotate with delim by ','.\nPossible list:"+','.join(possible_list)+'. (Example: "refGene,cytoBand,1000g2015aug_all,exac03,avsnp150") (Default=refGene)')
parser.add_argument("-b","--build", default='hg19', help="Build version of vcf (Default=hg19) If you choose mm10, only refGene will be annotated", choices=["hg19","hg38","mm10"])
parser.add_argument("-k","--keep", action="store_true", help="Don't remove intermediate files")
parser.add_argument("-s","--short_dbnsfp", action="store_true", help="Give summary version of dbnsfp when user want to dbnsfp annotation. Recommand to use this option because dbnsfp gives you too many columns")
args=parser.parse_args()

db_path='/home/users/sypark/03_Tools/annovar/humandb/'
if args.build == 'mm10':
	args.item='refGene'
	db_path='/home/users/sypark/03_Tools/annovar/mm10/'
	print('###Only refGene can be annotated for mm10###')


if 'dbnsfp35a' in args.items:
	item_list=args.items.split(',')
	del item_list[item_list.index('dbnsfp35a')]
	item_list.append('dbnsfp35a')
	args.items=','.join(item_list)

fn=args.vcf
input_list=args.items.split(',')
operation_list=[]
for name in input_list:
	if name not in possible_list:
		print("The input is not in possible list")
		print(name)
		sys.exit()
	if name == 'refGene':
		operation_list.append('g')
	elif name=='cytoBand':
		operation_list.append('r')
	
	else:
		operation_list.append('f')

os.system('/home/users/sypark/03_Tools/annovar/convert2annovar.pl --includeinfo -format vcf4old '+fn+' > '+fn+'.avinput')
os.system('/home/users/sypark/03_Tools/annovar/table_annovar.pl '+fn+'.avinput '+db_path+' -buildver '+args.build+' -out '+fn+'.anv -protocol '+args.items+' -operation '+','.join(operation_list)+' -remove -nastring . --dot2underline --otherinfo')

out_file=open(fn+'.anv','w')
in_file=open(fn)
in_line=in_file.readline().strip()
col_name1=''
while in_line[0]=='#':
	if in_line[0:4]=='#CHR':
		col_name1=in_line
	else:
		out_file.write(in_line+'\n')
	in_line=in_file.readline().strip()
in_file.close()

in_file=open(fn+'.anv.'+args.build+'_multianno.txt')
in_line=in_file.readline().strip()
in_indi=in_line.split('\t')
if 'dbnsfp' in args.items and args.short_dbnsfp == True:
	a=-72
	target_col=[a+3,a+6,a+9,a+12,a+15,a+18,a+21]
	col_name2='\t'.join(in_indi[5:-71])+'\t'+';'.join(in_indi[i] for i in target_col)
else:
	col_name2='\t'.join(in_indi[5:-1])
scol=len(in_indi)-1
if col_name1 != '':
	col_name=col_name1+'\t'+col_name2
	out_file.write(col_name+'\n')
in_line=in_file.readline().strip()
while in_line:
	in_indi=in_line.split('\t')
	if 'dbnsfp' in args.items and args.short_dbnsfp == True:
		a=scol-71
		target_col=[a+3,a+6,a+9,a+12,a+15,a+18,a+21]
		out_file.write('\t'.join(in_indi[scol:])+'\t'+'\t'.join(in_indi[5:a+1])+'\t'+';'.join(list(in_indi[i] for i in target_col))+'\n')
	else:
		out_file.write('\t'.join(in_indi[scol:])+'\t'+'\t'.join(in_indi[5:scol])+'\n')
	in_line=in_file.readline().strip()
in_file.close()

if args.keep:
	'blank'
else:
	os.system('rm '+fn+'.avinput')
	os.system('rm '+fn+'.anv.'+args.build+'_multianno.txt')

