# Make pipeline for github

stop('abc')

#load libraries
library(tidyverse)
library(readxl)

# extract examples for github
id_dt <- read_excel('/home/users/sypark/00_Project/06_LineageTracing/meta_data/example_id_for_github.xlsx')
early_mt_path <- '/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/perLineage/DB10_early_pointmts.txt'
el_dt <- read_tsv(early_mt_path, col_types = cols(`#CHROM`='c'))
early_ids <- el_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) %>% pull(var_id)

vs_path='/home/users/sypark/00_Project/06_LineageTracing/db10/03_varscan2'
hc_path='/home/users/sypark/00_Project/06_LineageTracing/db10/04_haplotypecaller'
save_path = '/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter/examples'

#Varscan sampling
for (i in 1:10){
  t_id <- id_dt$sample_id[i]
  ex_id <- id_dt$example_id[i]
  print(t_id)
  subs_path='/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution'
  passed_dt <- read_tsv(paste0(subs_path, '/',t_id,'.VS_HC_common.cmn_fi.readinfo.readc.rasm.snuN30.bgiN24.anv.fi'), col_types=cols(`#CHROM`='c'))
  passed5pct_ids <- passed_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) %>% sample_n(nrow(passed_dt)%/%20) %>%
    pull(var_id)
  dt=read_tsv(paste0(vs_path,'/',t_id,'.vs.snp.vcf.gz'), comment = '##', col_types = cols(`#CHROM`='c'))
  mdt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
  union_ids <- union(early_ids, passed5pct_ids)
  f1_dt <- mdt %>% filter(var_id %in% union_ids)
  nrow(f1_dt)
  f2_dt <- mdt %>% filter(!(var_id %in% union_ids)) %>% sample_n(10000 - nrow(f1_dt))
  sampled_dt <- bind_rows(f1_dt, f2_dt) %>% arrange(`#CHROM`, POS)
  write_tsv(sampled_dt, paste0(save_path, '/', ex_id,'.vs.snp.vcf'))
}

#Haplotype Caller sampling
for (i in 1:10){
  t_id <- id_dt$sample_id[i]
  ex_id <- id_dt$example_id[i]
  print(t_id)
  subs_path='/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution'
  passed_dt <- read_tsv(paste0(subs_path, '/',t_id,'.VS_HC_common.cmn_fi.readinfo.readc.rasm.snuN30.bgiN24.anv.fi'), col_types=cols(`#CHROM`='c'))
  passed5pct_ids <- passed_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) %>% sample_n(nrow(passed_dt)%/%20) %>%
    pull(var_id)
  dt=read_tsv(paste0(hc_path,'/',t_id,'.hc.snv.vcf'), comment = '##', col_types = cols(`#CHROM`='c'))
  mdt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
  union_ids <- union(early_ids, passed5pct_ids)
  f1_dt <- mdt %>% filter(var_id %in% union_ids)
  nrow(f1_dt)
  f2_dt <- mdt %>% filter(!(var_id %in% union_ids)) %>% sample_n(10000 - nrow(f1_dt))
  sampled_dt <- bind_rows(f1_dt, f2_dt) %>% arrange(`#CHROM`, POS)
  write_tsv(sampled_dt, paste0(save_path, '/', ex_id,'.hc.snv.vcf'))
}




#make union set of SNVs from 10 samples to make sampled panel of normals
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter', pattern='.cmn_fi$', full.names = T)
merged_dt <- do.call(rbind, lapply(file_list, function(x) read_tsv(x) %>% select(`#CHROM`, POS, REF, ALT)))
merged_dt <- unique(merged_dt)
merged_dt <- merged_dt %>% arrange(`#CHROM`, POS) %>% mutate(t_bin = POS%/%1000000, chr_tbin_id = paste0(`#CHROM`,'.',t_bin))
merged_dt %>% group_by(chr_tbin_id) %>% count() %>% nrow()

#extract line from snuN30 dataset
fd_path='/home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/01_SNU30/splitted_1Mbp/'
save_path='/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter/examples/panel_of_normal/'
library(doMC)
registerDoMC(8)
merged_ref <- foreach (t_id = merged_dt$chr_tbin_id %>% unique(), .combine = bind_rows) %dopar% {
  print(t_id)
  t_chr <- unlist(strsplit(t_id,'\\.'))[1]
  ref_dt <- read_tsv(paste0(fd_path,'N30_chr',t_chr,'.Q0q0.mpileup.call.calc.',t_id), 
                     col_names = c('#CHROM','POS','REF','Total','Abase','Gbase','Cbase','Tbase'),
                     col_types = cols(`#CHROM`='c'))
  target_pos <- merged_dt %>% filter(chr_tbin_id == t_id) %>% pull(POS)
  ref_dt <- ref_dt %>% filter(POS %in% target_pos) %>%mutate_at(vars(ends_with('base')), list(~ map_chr(., ~ paste(unlist(strsplit(.x,';'))[1:2], collapse=';'))))
  ref_dt
}
nrow(merged_ref)
merged_ref %>% write_tsv(paste0(save_path, 'NP1.txt'))

#extract line from bgiN24 dataset
fd_path = '/home/users/team_projects/LungAdeno_WGS_sypark/03_Normal_Panels/01_PointMts/04_BGI24/splitted_1Mbp/'
save_path='/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter/examples/panel_of_normal/'
library(doMC)
registerDoMC(8)
merged_ref <- foreach (t_id = merged_dt$chr_tbin_id %>% unique(), .combine = bind_rows) %dopar% {
  print(t_id)
  t_chr <- unlist(strsplit(t_id,'\\.'))[1]
  ref_dt <- read_tsv(paste0(fd_path,'bgiN24.24s.q0Q0.chr',t_chr,'.mpileup.snv.edit.',t_id,'.gz'), 
                     col_names = c('#CHROM','POS','REF','Total','Abase','Gbase','Cbase','Tbase'),
                     col_types = cols(`#CHROM`='c'))
  target_pos <- merged_dt %>% filter(chr_tbin_id == t_id) %>% pull(POS)
  ref_dt <- ref_dt %>% filter(POS %in% target_pos) %>%mutate_at(vars(ends_with('base')), list(~ map_chr(., ~ paste(unlist(strsplit(.x,';'))[1:2], collapse=';'))))
  ref_dt
}
nrow(merged_ref)
merged_ref %>% write_tsv(paste0(save_path, 'NP2.txt'))


#make bed for slicing BAMs

file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter', pattern='.cmn_fi$', full.names = T)
save_path='/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter/examples/'
for (i in 1:length(file_list)){
  print(file_list[i])
  r <- regexec(".*/(.*).VS",file_list[i])
  m <- regmatches(file_list[i], r)
  sampleid <- map_chr(m,2)
  read_tsv(file_list[i], col_types = cols(`#CHROM`='c')) %>% 
    mutate(start_pos = POS-1) %>% select(`#CHROM`, start_pos, POS) %>% 
    write_tsv(paste0(save_path, sampleid,'.bed'), col_names = F)
}


#make union bed for control BAM slicing
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter/examples', pattern = '.bed$', full.names = T)
length(file_list)
read_tsv(file_list[1])
dt <- do.call(rbind, lapply(file_list, function(x) read_tsv(x, col_names=c('chrom','start','end'), col_types = cols(chr='c'))))
dt$chrom <- factor(dt$chrom, levels = c(1:22,'X'))
dt %>% arrange(chrom, start) %>% unique() %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/git_hub/Substitution_filter/examples/example1to10union.bed', col_names = F)

