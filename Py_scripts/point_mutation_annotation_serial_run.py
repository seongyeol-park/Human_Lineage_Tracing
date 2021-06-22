#191005 written

import argparse, sys, subprocess
parser = argparse.ArgumentParser(
	description = "this script run serial annotation scripts for point mutations. 1: readinfo, 2: read count of paired normal, 3: read local reassembly info")
parser.add_argument("-i","--input", type=str, help="input vcf file")
parser.add_argument("-t","--tumor_bam", type=str, help="input tumor bam file")
parser.add_argument("-n","--normal_bam", type=str, help="input normal bam file")
parser.add_argument("-r","--run", type=int, help= "starting number of process which you would like to run. 1: 1-3 scripts, 2: 2-3 scripts only, 3: 3 script only. (default = 1)", default = 1 )
args=parser.parse_args()
try: 
	print(','.join([args.input, args.tumor_bam, args.normal_bam, str(args.run)]))
except:
	print('ERROR: false arguments')
	parser.print_help(sys.stderr)
	sys.exit(1)


scrDir='/home/users/team_projects/Point_Mutation_Annotation_Toolkit'
scr1='readinfo_anno_bwa_190920.py'
scr2='readcount_only_anno_bwa_190729.py'
scr3='read_local_reassembly_190729.py'

cmd1=' '.join(['python', scrDir+'/'+scr1, args.input, args.tumor_bam])
cmd2=' '.join(['python', scrDir+'/'+scr2, args.input+'.readinfo', args.normal_bam])
cmd3=' '.join(['python', scrDir+'/'+scr3, args.input+'.readinfo.readc', args.tumor_bam, args.normal_bam])
rm1='rm -f '+args.input+'.readinfo'
rm2='rm -f '+args.input+'.readinfo.readc'

cmd_list=[cmd1, cmd2, cmd3, rm1, rm2]
merged_cmd=' && '.join(cmd_list[args.run-1:])
print(merged_cmd)
res=subprocess.Popen(merged_cmd, shell=True)
if res.wait() != 0:
	output, error = res.communicate()
	print(error)
	sys.exit(1)
