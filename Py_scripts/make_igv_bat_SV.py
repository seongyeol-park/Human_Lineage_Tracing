#191015 exit when there is no SV, goto chr1 before making start.png
#191030 write samplie id at file name rather than folder name
#200420 make igv_bat directory and save output files into the folder

import sys,os, subprocess
dist=1500
print(sys.argv[1])
in_file=open(sys.argv[1])   # SV list
t_bam=sys.argv[2]  #absolute path
n_bam=sys.argv[3]  #absolute path
gff=sys.argv[4] #absolute path
save_dir=sys.argv[5]  #absolute path
sampleid=sys.argv[6]

print('recommendation: set insert size option 50, 1000')

if save_dir[-1]=='/':
	save_dir=save_dir[0:-1]

os.system('mkdir -p igv_bat')


#check line_number
res=subprocess.Popen('cat '+sys.argv[1]+' | grep -v ^# | wc -l', shell = True, stdout=subprocess.PIPE)
res.wait()
line_num = int(res.stdout.readline().rstrip())
print(line_num)
if line_num == 0:
	print('no SV was found')
	sys.exit()

os.system('mkdir -p '+save_dir)
out_file=open('igv_bat/'+sys.argv[1]+'.igv.bat','w')
out_file.write('new\n')
out_file.write('snapshotDirectory '+save_dir+'\n') # absolute path
out_file.write('load '+t_bam+'\n') #absolute path
out_file.write('load '+n_bam+'\n') #absolute path
out_file.write('load '+gff+'\n') #absolute path
out_file.write('collapse\n')
gff_track=gff.split('/')[-1]
out_file.write('expand '+gff_track+'\n')
out_file.write('setSleepInterval 500ms\n')
out_file.write('goto chr1\n')
out_file.write('snapshot start.png\n')

#make bat
in_line=in_file.readline().strip()
while in_line:
	if in_line[0]=='#':
		'blank'
	else:
		in_indi=in_line.split('\t')
		chr1=in_indi[0]
		pos1=in_indi[1]
		chr2=in_indi[2]
		pos2=in_indi[3]
		terinfo=in_indi[5]
		svtype=in_indi[6]
		precise=in_indi[7]
		upper=int(pos1)-dist
		lower=int(pos1)+dist
		out_file.write('goto chr'+chr1+':'+str(upper)+'-'+str(lower)+'\n')
		out_file.write('snapshot '+sampleid+'.'+'_'.join([chr1,pos1,chr2,pos2])+'_'+terinfo+'_'+svtype+'_'+precise+'_1.png\n')
		if chr2 != '.':
			upper=int(pos2)-dist
			lower=int(pos2)+dist
			out_file.write('goto chr'+chr2+':'+str(upper)+'-'+str(lower)+'\n')
			out_file.write('snapshot '+sampleid+'.'+'_'.join([chr1,pos1,chr2,pos2])+'_'+terinfo+'_'+svtype+'_'+precise+'_2.png\n')
	in_line=in_file.readline().strip()

