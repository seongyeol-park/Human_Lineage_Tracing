#probability of mutation overlap by chance simple statistics #reviewer2 number2

library(tidyverse)
library(ggsci)
library(cowplot)

meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210330.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210227.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_200918.txt')

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

source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Lineage_Tracing_Color_Palette.R')



Assign_Lineage_ID <- function(ts_dt, this_db){
  merged_tbl = tibble()
  file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/'), pattern = 'L.*\\.pointmts.txt$', full.names = T)
  for (n in 1:length(file_list)){
    dt <- read_tsv(file_list[n], col_types = cols(`#CHROM`='c', REF='c', ALT='c'))
    l_id <- gsub('.pointmts.txt','',unlist(strsplit(file_list[n],'//'))[2])
    dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'), lineage_id2 = paste(this_db, l_id, sep='_')) %>%
      separate("blood_dp;blood_var;blood_vafpct", c("blood_dp","blood_var","blood_vafpct"), sep=';') %>%
      mutate(WGS_bloodVAF = as.numeric(blood_vafpct)/100) %>% select(var_id, lineage_id2, WGS_bloodVAF)
    merged_tbl <- bind_rows(merged_tbl, left_join(ts_dt, dt) %>% filter(is.na(lineage_id2) == F))
  }
  return(merged_tbl)
}


#find overlapped mutations and take bamsnap
blood_tbl=tribble(~deadbody, ~blood_bam,
                  "DB2","2_blood_3s_merged.s.md.ir.br.bam",
                  "DB3","3_blood_4s_merged.bam",
                  "DB5","5_blood_3s_merged.s.md.ir.br.bam",
                  "DB6","6_Blood-merged.rg.bam",
                  "DB8","8_Blood_3s_merged.s.md.ir.br.bam",
                  "DB9","9_Blood_3s_merged.s.md.ir.br.bam",
                  "DB10","10_Blood_3s_merged.s.md.ir.br.bam")
                  

if(F){  #you don't need to do this again
  #make table and bamsnap for visual insepction of mutaitons assigned to multiple lineages
  if(F){
    merged_dt <- tibble()
    for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
      print(this_db)
      fd_path = '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/'
      fn= paste0(this_db, '_pointmt_merged.txt')
      dt <- read_tsv(paste0(fd_path, fn), col_types = cols(`#CHROM`='c'))
      dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) 
      f_dt <- dt %>% filter(caseNo> 1) 
      mf_dt <- Assign_Lineage_ID(f_dt, this_db)
      mf_dt <- left_join(mf_dt, lng_dt %>% select(lineage_id2, time_of_lineage, branch_n))
      #mf_dt %>% group_by(var_id) %>% dplyr::count() %>% filter(n>1)
      overlapped_ids <- mf_dt %>% group_by(var_id) %>% dplyr::count() %>% filter(n>1) %>% pull(var_id)
      fmf_dt <- mf_dt %>% filter(var_id %in% overlapped_ids)
      fmf_dt <- fmf_dt %>% mutate(deadbody = this_db, lineage_id = lapply(lineage_id2, function(x) unlist(strsplit(x, '_'))[2]) %>% as.character(),
                                  var_id2 = gsub('_',":",var_id))
      merged_dt <- bind_rows(merged_dt, fmf_dt %>% group_by(deadbody, var_id2, lineage_id) %>% count() %>% select(-n))
      
      if(F){ #bamsanp
        for (i in 1:nrow(fmf_dt)){
          print(i)
          bamsnap_path= '/home/users/sypark/anaconda3/bin/bamsnap'
          this_chrom = fmf_dt$`#CHROM`[i]
          this_pos = fmf_dt$POS[i]
          this_caselist = fmf_dt$case_list[i]
          this_varid = fmf_dt$var_id[i]
          bam_fd = paste0('/home/users/team_projects/Lineage_tracing/',this_db,'/02_BAM/')
          bams <- paste0(bam_fd, unlist(strsplit(this_caselist,';')), '.s.md.ir.br.bam')
          blood_bam = paste0(bam_fd,blood_tbl$blood_bam[blood_tbl$deadbody == this_db])
          bams <- c(bams, blood_bam)
          
          l_list=c()
          for (t_id in unlist(strsplit(this_caselist,';'))){
            print(t_id)
            l_list = c(l_list, meta_dt$lineage_id[meta_dt$sample_id == t_id])
          }
          l_list <- c(l_list, 'blood')
          l_list <- paste0(l_list, ',',this_varid)
          
          snap_pos = paste0(this_chrom,':',this_pos)
          out_fd = '/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/bamsnap_overlapped_mutations'
          out_file= paste0(this_db, '_', this_varid,'.png')
          cmd = paste0(bamsnap_path,' -bam ', paste(bams, collapse=" "),' -refversion hg19 -pos ',snap_pos,' -out ',out_fd, '/',out_file, ' -title ',paste(l_list, collapse=' '),' -show_soft_clipped -margin 200')
          print(cmd)
          system(cmd)
        }
      }
    }
    merged_dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/overlapped_mutations.tsv')  
  }
  
  #make perSample link after regeneration of mutation list
  if(F){
    for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
      print(this_db)
      id_list <- meta_dt %>% filter(deadbody == this_db) %>% pull(sample_id)
      fd_path = paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perSample/')
      out_fd_path = '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/'
      cmd_list <- paste0('ln -s ',fd_path, id_list,'.edited.vcf.anv ',out_fd_path)
      lapply(cmd_list, system)
    }
  }
  
}



#count mutations assinged to multiple lineages after visual inspection and regeneration of the list
if(F){
  total_tbl <- tibble(deadbody = as.character(), n_unique = as.numeric(), n_shared_by_late_br = as.numeric(), n_shared_early = as.numeric(), n_shared_recurrent = as.numeric())
  total_recur <- tibble()
  total_per_sample <- tibble()
  this_db = 'DB6'
  n=0
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n=n+1
    print(this_db)
    id_list <- meta_dt %>% filter(deadbody == this_db) %>% pull(sample_id)
    
    i=1
    merged_dt <- tibble()
    for (i in 1:length(id_list)){
      print(i)
      this_id=id_list[i]
      file_path = paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody==this_db],'/perSample/', this_id, '.pointmt')
      dt <- read_tsv(file_path, col_types = cols(`#CHROM`='c'))
      #colnames(dt) <- c('#CHROM','POS','ID','REF','ALT','info')
      #dt <- dt %>% separate(info, c('ref','var','ukn','vafpct'), sep=';') %>% mutate(sample_id = this_id)
      dt <- dt %>% mutate(sample_id = this_id)
      merged_dt <- bind_rows(merged_dt, dt)
    }
    merged_dt <- merged_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
    nrow(merged_dt)
    m_dt <- Assign_Lineage_ID(merged_dt, this_db)
    nrow(m_dt)
    m_dt <- left_join(m_dt, lng_dt %>% select(lineage_id2, time_of_lineage, branch_n))
    tmp1 <- m_dt %>% select(sample_id, var_id) %>% unique() %>% group_by(var_id) %>% summarise(sample_n = n())
    tmp2 <- m_dt %>% select(var_id, lineage_id2) %>%  unique() %>% group_by(var_id) %>% summarise(lineage_n=n())
    m_dt <- left_join(left_join(m_dt, tmp1),tmp2)
    
    tmp0 <- merged_dt %>% group_by(sample_id) %>% summarise(n_total = n())
    tmp1 <- m_dt %>% filter(sample_n == 1) %>% group_by(sample_id) %>% summarise(n_unique = n())
    tmp2 <- m_dt %>% filter(sample_n > 1 & lineage_n == 1) %>% filter(time_of_lineage != 'early') %>% group_by(sample_id) %>% summarise(n_shared_by_late_br = n())
    tmp3 <- m_dt %>% filter(sample_n > 1 & lineage_n == 1 & time_of_lineage == 'early') %>% group_by(sample_id) %>% summarise(n_shared_early = n())
    tmp4 <- m_dt %>% filter(sample_n >1 & lineage_n >1)%>% select(var_id, sample_id) %>% unique() %>% group_by(sample_id) %>% summarise(n_shared_recur = n())
    tmp <- left_join(left_join(left_join(left_join(tmp0, tmp1), tmp2), tmp3), tmp4)
    tmp[is.na(tmp)] <- 0
    tmp$deadbody = this_db
    total_per_sample <- bind_rows(total_per_sample, tmp)
    
    n_unique <- m_dt %>% filter(sample_n == 1) %>% pull(var_id) %>% unique() %>% length()
    n_shared_by_late_br <- m_dt %>% filter(sample_n > 1 & lineage_n == 1) %>% filter(time_of_lineage != 'early') %>% pull(var_id) %>% unique() %>% length()
    n_shared_early <- m_dt %>% filter(sample_n > 1 & lineage_n == 1 & time_of_lineage == 'early') %>% pull(var_id) %>% unique() %>% length()
    n_shared_recur <- m_dt %>% filter(sample_n >1 & lineage_n >1)%>% select(var_id, lineage_id2) %>% unique() %>% nrow()
    m_dt %>% filter(sample_n >1 & lineage_n >1) %>% write_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/', this_db, '.recurret_mutations.tsv'))
    total_tbl[n,] <- list(this_db, n_unique, n_shared_by_late_br, n_shared_early, n_shared_recur)
    total_recur <- bind_rows(total_recur, m_dt %>% filter(sample_n >1 & lineage_n >1))
  }
  total_per_sample %>% mutate(calc = n_unique + n_shared_by_late_br + n_shared_early + n_shared_recur) %>% mutate(gap = n_total - calc) %>% pull(gap) %>% summary()
  total_per_sample %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/pointmt_count_unique_shared_perSample.tsv')
  total_tbl <- total_tbl %>% mutate(n_total = n_unique + n_shared_by_late_br + n_shared_early + n_shared_recurrent)
  total_tbl <- total_tbl %>% mutate(n_nonEarly = n_total - n_shared_early)
  total_tbl %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/pointmt_count_shared_nonshared.tsv')
  
}
total_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/pointmt_count_shared_nonshared.tsv')
sum(total_tbl$n_unique)/sum(total_tbl$n_nonEarly)
sum(total_tbl$n_shared_by_late_br)/sum(total_tbl$n_nonEarly)
sum(total_tbl$n_shared_recurrent)/sum(total_tbl$n_nonEarly)
sum(total_tbl$n_shared_recurrent)
total_tbl <- left_join(total_tbl, DB_dt %>% select(deadbody, sample_n))
pdf("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig10a.pdf")
ggplot(total_tbl, aes(x=n_total, y=n_shared_recurrent))+
  geom_point(size=5, aes(color=deadbody), alpha=0.8)+
  scale_color_manual(values=db_pal, name = "")+
  xlab("No. of total mutations")+ylab("No. of recurrent mutations")+
  theme_syp
dev.off()
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/', pattern = 'recurret_mutations.tsv$', full.names = T)
merged_recur <- do.call(rbind, lapply(file_list, function(x) read_tsv(x)))
merged_recur$var_id %>% unique() %>% length()
merged_recur %>% select(`#CHROM`, POS, ID, REF, ALT) %>% unique() %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/unique_set_recurrent_mutations.tsv')


#reload
total_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/pointmt_count_shared_nonshared.tsv')

total_late <- sum(total_tbl$n_total - total_tbl$n_shared_early)
total_late
unique_sum <- total_tbl$n_unique %>% sum()
unique_sum/total_late
shared_late_sum <- total_tbl$n_shared_by_late_br %>% sum()
shared_late_sum/total_late
shared_recur <- total_tbl$n_shared_recurrent %>% sum()
shared_recur/total_late

DB6_total <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/SBS15718_merged/DB6_pointmt_merged.txt.signature_exposures.SBS15718.tsv')
DB6_early <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/perLineage/sig_subs_SBS15718_a1/DB6_early_SNVs.txt.signature_exposures.tsv')
DB6_recurrent <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/SBS15718/DB6.recurret_mutations.tsv.signature_exposures.SBS15718.tsv')

DB10_total <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/SBS15718_merged/DB10_pointmt_merged.txt.signature_exposures.SBS15718.tsv')
DB10_early <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/perLineage/sig_subs_SBS15718_a1/DB10_early_SNVs.txt.signature_exposures.tsv')
DB10_recurrent  <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/n_recurrent_mutations/SBS15718/DB10.recurret_mutations.tsv.signature_exposures.SBS15718.tsv')



DB6_total <- DB6_total %>% dplyr::rename(signature = `Signature #`) %>% filter(signature != 'Unexplained')
DB6_total <- DB6_total %>% mutate(DB6_total_prop = Exposure*100/sum(Exposure)) %>% select(-Exposure, -Proportion)

DB6_early <- DB6_early %>% dplyr::rename(signature = `Signature #`) %>% filter(signature != 'Unexplained')
DB6_early <- DB6_early %>% mutate(DB6_early_prop = Exposure*100/sum(Exposure)) %>% select(-Exposure, -Proportion)

DB6_recurrent <- DB6_recurrent %>% dplyr::rename(signature = `Signature #`) %>% filter(signature != 'Unexplained')
DB6_recurrent <- DB6_recurrent %>% mutate(DB6_recurrent_prop = Exposure*100/sum(Exposure))%>% select(-Exposure, -Proportion)

DB10_total <- DB10_total %>% dplyr::rename(signature = `Signature #`) %>% filter(signature != 'Unexplained')
DB10_total <- DB10_total %>% mutate(DB10_total_prop = Exposure*100/sum(Exposure)) %>% select(-Exposure, -Proportion)

DB10_early <- DB10_early %>% dplyr::rename(signature = `Signature #`) %>% filter(signature != 'Unexplained')
DB10_early <- DB10_early %>% mutate(DB10_early_prop = Exposure*100/sum(Exposure)) %>% select(-Exposure, -Proportion)

DB10_recurrent <- DB10_recurrent %>% dplyr::rename(signature = `Signature #`) %>% filter(signature != 'Unexplained')
DB10_recurrent <- DB10_recurrent %>% mutate(DB10_recurrent_prop = Exposure*100/sum(Exposure)) %>% select(-Exposure, -Proportion)


DB6_sig_dt <- left_join(left_join(DB6_total, DB6_early), DB6_recurrent)
DB6_sig_dt <- DB6_sig_dt %>% gather(-signature, key="type",value="proportion")
DB6_sig_dt$type <- factor(DB6_sig_dt$type, levels = c('DB6_total_prop','DB6_early_prop','DB6_recurrent_prop'))
g1 <- ggplot(DB6_sig_dt, aes(x=type))+
  geom_bar(aes(y=proportion, fill=signature), stat="identity")+
  scale_fill_manual(values=sig_pal)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels = c('Total','Early','Recurrent'))+
  theme_syp + theme(axis.title.x = element_blank(), legend.position='none')+
  ggtitle('DB6')

DB10_sig_dt <- left_join(left_join(DB10_total, DB10_early), DB10_recurrent)
DB10_sig_dt <- DB10_sig_dt %>% gather(-signature, key="type",value="proportion")
DB10_sig_dt$type <- factor(DB10_sig_dt$type, levels = c('DB10_total_prop','DB10_early_prop','DB10_recurrent_prop'))
g2 <- ggplot(DB10_sig_dt, aes(x=type))+
  geom_bar(aes(y=proportion, fill=signature), stat="identity")+
  scale_fill_manual(values=sig_pal)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(labels = c('Total','Early','Recurrent'))+
  theme_syp + theme(axis.title.x = element_blank())+
  ggtitle('DB10')

plot_grid(g1,g2, nrow=1, align='h', axis='tb', rel_widths = c(1,1.2))


#recount recurrent muatations

this_db='DB6'
lng_dt %>% filter(is.na(n_pointmt) == F) %>%  filter(time_of_lineage != 'early') %>% pull(n_pointmt)%>% sum()  #total non early mutations
lng_dt %>% filter(is.na(n_pointmt) == F & deadbody == this_db) %>%  filter(time_of_lineage == 'early') 
lng_dt %>% filter(is.na(n_pointmt) == F & deadbody == this_db) %>%  filter(time_of_lineage != 'early') %>% pull(n_pointmt)%>% sum()  #total non early in this_db
lng_dt %>% filter(is.na(n_pointmt) == F) %>%  filter(time_of_lineage != 'early') 
