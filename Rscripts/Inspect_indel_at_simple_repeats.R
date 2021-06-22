#Find MSI related indels
stop('abc')

library(tidyverse)


pointmt_path_tbl <- tribble(
  ~deadbody, ~pointmt_path,
  "DB2","/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/",
  "DB3","/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/",
  "DB5","/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/",
  "DB6","/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/",
  "DB8","/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/",
  "DB9","/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/",
  "DB10","/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/"
)

lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210330.txt') %>% filter(current_final_for_lineage == 'Y')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_210322.txt')

#inspect whether there is a sample with high indels burden
m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with('SBS')), list(~ n_snv*.))
m_meta_dt <- m_meta_dt %>% mutate(SBS7 = SBS7a + SBS7b + SBS7c + SBS7d) %>% filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10'))
ggplot(m_meta_dt, aes(x=SBS7, y=n_indel))+
  geom_point()+
  facet_wrap(~deadbody)



#merge indel results and find indels occured at microsatellite region
merged_dt <- tibble()
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  print(this_db)
  fd_path = paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db], 'perLineage/')
  file_name = paste0(this_db,'_early_indels.txt.rpct')
  file_name2 = paste0(this_db, '_lateNonUV_indels.txt.rpct')
  file_name3 = paste0(this_db, '_lateUV_indels.txt.rpct')
  
  indel_dt1 <- read_tsv(paste0(fd_path, file_name), col_types= cols(`#CHROM`='c', ALT='c'))
  indel_dt2 <- read_tsv(paste0(fd_path, file_name2), col_types= cols(`#CHROM`='c',ALT='c'))
  indel_dt3 <- read_tsv(paste0(fd_path, file_name3), col_types= cols(`#CHROM`='c',ALT='c'))
  
  indel_dt <- bind_rows(bind_rows(indel_dt1, indel_dt2), indel_dt3)
  indel_dt <- indel_dt %>% separate("repeat_unit;ref_repeat_count", c("repeat_unit","ref_repeat_count"), sep=';', convert=T) %>%
    mutate(repeat_size= nchar(repeat_unit))
  merged_dt <- bind_rows(merged_dt, indel_dt)
}

merged_dt <- left_join(merged_dt, lng_dt)

merged_dt %>% filter(repeat_size < 7 & ref_repeat_count >=5) %>% group_by(deadbody) %>% dplyr::count()
merged_dt %>% filter(repeat_size < 7 & ref_repeat_count >=5) %>%
  filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10')) %>% nrow()

merged_dt %>% filter(repeat_size < 7 & ref_repeat_count >=5 & time_of_lineage == 'early') %>% group_by(deadbody) %>% dplyr::count()
merged_dt %>% filter(repeat_size < 7 & ref_repeat_count >=5 & time_of_lineage == 'early') %>%
  filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10')) %>% nrow()




f_dt <- merged_dt %>% filter(repeat_size < 7 & ref_repeat_count >=5 & time_of_lineage == 'early') %>%
  filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10'))

f_dt %>% select(1:11) %>% arrange(lineage_id2) %>% View()

tmp_dt1 <- lng_dt %>% filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10') & time_of_lineage == 'early') %>% group_by(deadbody) %>% summarise(early_point_mutations = sum(n_pointmt), early_indels = sum(n_indels))
tmp_dt2 <- f_dt %>% group_by(deadbody) %>% summarise(early_indels_in_MS_region = n())
tmp_dt <- left_join(tmp_dt1, tmp_dt2)
tmp_dt[is.na(tmp_dt)] <- 0
tmp_dt <- left_join(tmp_dt, DB_dt %>% select(deadbody, sample_n))
tmp_dt <- tmp_dt %>% dplyr::rename(n_clonal_samples = sample_n, donor_id = deadbody) %>% select(donor_id, n_clonal_samples, everything())
tmp_dt %>% mutate(prop = early_indels_in_MS_region/early_indels)



#heatmap of DB3 indel
db3_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/perLineage/DB3_early_pointmts.txt.43sCall')
db3_indel_dt <- db3_dt %>% filter(nchar(REF) != 1 | nchar(ALT) != 1)
colnames(db3_indel_dt)
mx <- db3_indel_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) %>% select(var_id, ends_with('vafpct'), -starts_with('blood')) %>% as.data.frame() %>%
  column_to_rownames('var_id') %>% as.matrix()
library(pheatmap)
pheatmap(mx)



#bamsnap
bamsnap_path= '/home/users/sypark/anaconda3/bin/bamsnap'
this_chrom = '16'
this_pos = 55673730
snap_pos = paste0(this_chrom,':',this_pos)
out_fd =  "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/Inspect_indel_at_simple_repeat/"
bam_fd = '/home/users/team_projects/Lineage_tracing/DB3/02_BAM/'

#L1 samples
out_file= paste0(this_chrom, '_',this_pos, '.L1.bamsnap.png')
sample_ids = c('3_RLM_001E10','3_RCS_002A2')
bam_path = paste0(bam_fd, sample_ids, '.s.md.ir.br.bam')
titles <-  c('3_RLM_001E10','3_RCS_002A2')
#L2 samples
out_file= paste0(this_chrom, '_',this_pos, '.L2.bamsnap.png')
sample_ids = c('3_RT_Fb_001B1','3_LLL-1_Fb_001C3')
bam_path = paste0(bam_fd, sample_ids, '.s.md.ir.br.bam')
titles <- c('3_RT_Fb_001B1','3_LLL-1_Fb_001C3')

cmd = paste0(bamsnap_path,' -bam ', paste(bam_path, collapse=' '), ' -pos ',snap_pos,' -out ',out_fd, '/',out_file, ' -title ',paste(titles, collapse=' '),' -margin 50 -refversion hg19')
print(cmd)
system(cmd)

img <- readPNG(paste0(out_fd, '/',out_file))
plot.new()
grid::grid.raster(img)
