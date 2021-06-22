#!/bin/bash

#2021-03-15 change path of snuN30

input=$1

echo $input && echo "snuN30" && (python /home/users/sypark/01_Python_files/normal_matrix/AddNpanelToVCFsnp_vaf.py $input /home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/01_SNU30/snuN30.30s.q0Q0.chr1.mpileup.snv.edit.gz snuN30) || { c=$?;exit $c; } && echo "done" && echo "bgiN24" &&
(python /home/users/sypark/01_Python_files/normal_matrix/AddNpanelToVCFsnp_vaf.py $input.snuN30 /home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/04_BGI24/bgiN24.24s.q0Q0.chr1.mpileup.snv.edit.gz bgiN24) || { c=$?;exit $c; } && echo "done"
