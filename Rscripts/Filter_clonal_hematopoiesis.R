library(tidyverse)

#Filter blood specific variants
#previous filter
if(F){
  path_tbl <- tribble(
    ~deadbody, ~path, ~path2,
    "DB2", '/home/users/sypark/00_Project/06_LineageTracing/db2/11_clonal_hematopoiesis/HC_vcf_link/DB2_HC_SNV_bloodOnly.txt.readinfo.readc','/home/users/sypark/00_Project/06_LineageTracing/db2/11_clonal_hematopoiesis/HC_vcf_link/DB2_HC_SNV_bloodOnly.txt.readinfo.readc.fi.28sCall',
    "DB3", '/home/users/sypark/00_Project/06_LineageTracing/db3/15_clonal_hematopoiesis/HC_vcf_link/DB3_HC_SNV_bloodOnly.txt.readinfo.readc','/home/users/sypark/00_Project/06_LineageTracing/db3/15_clonal_hematopoiesis/HC_vcf_link/DB3_HC_SNV_bloodOnly.txt.readinfo.readc.fi.44sCall',
    "DB5",'/home/users/sypark/00_Project/06_LineageTracing/db5/12_clonal_hematopoiesis/HC_vcf_links/DB5_HC_SNV_bloodOnly.txt.readinfo.readc','/home/users/sypark/00_Project/06_LineageTracing/db5/12_clonal_hematopoiesis/HC_vcf_links/DB5_HC_SNV_bloodOnly.txt.readinfo.readc.fi.29sCall',
    "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/HC_vcf_links/DB6_HC_SNV_bloodOnly.txt.readinfo.readc','/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/HC_vcf_links/DB6_HC_SNV_bloodOnly.txt.readinfo.readc.fi.116sCall',
    "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/10_clonal_hematopoiesis/HC_vcf_link/DB8_HC_SNV_bloodOnly.txt.readinfo.readc','/home/users/sypark/00_Project/06_LineageTracing/db8/10_clonal_hematopoiesis/HC_vcf_link/DB8_HC_SNV_bloodOnly.txt.readinfo.readc.fi.48sCall',
    "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/10_clonal_hematopoiesis/HC_vcf_link/DB9_HC_SNV_bloodOnly.txt.readinfo.readc','/home/users/sypark/00_Project/06_LineageTracing/db9/10_clonal_hematopoiesis/HC_vcf_link/DB9_HC_SNV_bloodOnly.txt.readinfo.readc.fi.36sCall',
    "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/10_clonal_hematopoiesis/HC_vcf_link/DB10_HC_SNV_bloodOnly.txt.readinfo.readc','/home/users/sypark/00_Project/06_LineageTracing/db10/10_clonal_hematopoiesis/HC_vcf_link/DB10_HC_SNV_bloodOnly.txt.readinfo.readc.fi.40sCall'
  )
  
  false_list=list()
  false_list[["DB3"]] <- c('5_63662226_T_C')
  false_list[["DB5"]] <- c("3_174085929_G_A", "8_38467585_T_C","8_135916612_G_A", "12_115760614_C_T", "16_48766455_A_C", "16_48766457_A_C")
  false_list[["DB6"]] <- c("22_38083286_A_C", "X_132887787_A_C","X_132887786_C_G")
  false_list[["DB9"]] <- c("4_8541410_T_G")
  
  # filtering
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    print(this_db)
    dt <- read_tsv(path_tbl$path[path_tbl$deadbody == this_db],
                   col_types = cols(`#CHROM`='c'))
    colnames(dt)[(ncol(dt)-9):ncol(dt)] <- paste0('S',c(1:10))
    m_dt <- dt %>% separate(`ref_minMQ;med;max`, c("ref_minMQ","ref_medMQ","ref_maxMQ"), convert= T, sep =';') %>%
      separate(`var_minMQ;med;max`, c("var_minMQ","var_medMQ","var_maxMQ"), covert = T, sep =';') %>%
      separate(`ref_med;meanBQ`, c("ref_medBQ","ref_meanBQ"), convert = T, sep=';') %>%
      separate(`var_med;meanBQ`, c("var_medBQ","var_meanBQ"), convert =T, sep=';') %>%
      separate(`var_minMismatch;med;max`, c("var_minMismatch","var_medMismatch","var_maxMisMatch"), convert = T, sep=';') %>%
      separate(`refClip;pct;varClip;pct`, c("refClipN","refClipPct","varClipN","varClipPct"), convert = T, sep=';') %>%
      mutate_at(vars(starts_with("var_")), as.numeric) %>%
      separate(S1, c('S1_ref','S1_var','S1_vaf'),convert = T, sep=';') %>%
      separate(S2, c('S2_ref','S2_var','S2_vaf'),convert = T, sep=';') %>%
      separate(S3, c('S3_ref','S3_var','S3_vaf'),convert = T, sep=';') %>%
      separate(S4, c('S4_ref','S4_var','S4_vaf'),convert = T, sep=';') %>%
      separate(S5, c('S5_ref','S5_var','S5_vaf'),convert = T, sep=';') %>%
      separate(S6, c('S6_ref','S6_var','S6_vaf'),convert = T, sep=';') %>%
      separate(S7, c('S7_ref','S7_var','S7_vaf'),convert = T, sep=';') %>%
      separate(S8, c('S8_ref','S8_var','S8_vaf'),convert = T, sep=';') %>%
      separate(S9, c('S9_ref','S9_var','S9_vaf'),convert = T, sep=';') %>%
      separate(S10, c('S10_ref','S10_var','S10_vaf'),convert = T, sep=';') %>%
      mutate(VAF = (var_readN)/(var_readN + ref_readN)) %>%
      mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
    m_dt <- m_dt %>% mutate(NP_maxVAF = apply(m_dt[paste0('S',c(1:10),'_vaf')],1, function(x) max(x)))
    dim(m_dt)
    subdt <- m_dt %>% filter(var_medMQ >= 25 & 
                               ref_medMQ - var_medMQ < 10 & 
                               var_meanBQ >= 20 & 
                               NP_maxVAF < 1 & 
                               var_readN > 2 &
                               ref_minMQ > 0 & 
                               varClipPct < 95 &
                               !(var_id %in% false_list[[this_db]]) & 
                               var_medMismatch < 5)
    dim(subdt)
    write_tsv(subdt,paste0(path_tbl$path[path_tbl$deadbody == this_db],'.fi'))
  }
  
  #plot VAF heatmap for check
  library(ComplexHeatmap)
  for (this_db in c('DB2','DB5','DB6','DB8','DB9','DB10')){
    print(this_db)
    dt <- read_tsv(path_tbl$path2[path_tbl$deadbody == this_db],
                   col_types = cols(`#CHROM`='c'))
    mx <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) %>%
      select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
    print(Heatmap(mx))
  }
  
  # plot histogram
  p = list()
  n=0
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n = n+1
    print(this_db)
    dt <- read_tsv(paste0(path_tbl$path[path_tbl$deadbody == this_db], '.fi'),
                   col_types = cols(`#CHROM`='c'))
    p[[n]] <- ggplot(dt)+
      geom_histogram(aes(VAF), binwidth=0.02)+
      scale_x_continuous(limits = c(0,1), breaks=seq(0,1,0.1))+
      xlab("Variant-allele frequency")+
      ylab("No. of blood-specific subs")+
      ggtitle(paste0(this_db, ' (n=',nrow(dt),')'))
  }
  library(cowplot)
  plot_grid(plotlist = p, nrow=2)
}

#IGV view
subdt 
n=0

if(T){
  n=n+1
  print(n)
  print(paste(as.character(subdt[n,1:5])[c(1,2,4,5)], collapse='_'))
  system(paste('echo "goto', subdt[n,"#CHROM"], subdt[n,"POS"]-250,subdt[n,"POS"]+250,'" | nc localhost 60168', sep=' '))
  system('echo "sort" | nc localhost 60168')
  }

subdt[n,] %>% View()

##########filter again with more variants
#include some shared variants for validation
path_tbl2 <- tribble(
  ~deadbody, ~folder_path, ~file_name1,~blood_id,
  "DB2", '/home/users/sypark/00_Project/06_LineageTracing/db2/11_clonal_hematopoiesis/HC_vcf_link/','DB2_28s_HC_SNV_merged.txt.inc2_blood_3s_merged.less5.readinfo','2_blood_3s_merged',
  "DB3", '/home/users/sypark/00_Project/06_LineageTracing/db3/15_clonal_hematopoiesis/HC_vcf_link/','DB3_44s_HC_SNV_merged.txt.inc3_blood_4s_merged.less5.readinfo','3_blood_4s_merged',
  "DB5",'/home/users/sypark/00_Project/06_LineageTracing/db5/12_clonal_hematopoiesis/HC_vcf_links/','DB5_29s_HC_SNV_merged.txt.inc5_blood_3s_merged.less5.readinfo','5_blood_3s_merged',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/HC_vcf_links/','DB6_116s_HC_SNV_merged.txt.inc6_Blood-merged.less5.readinfo','6_Blood-merged',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/10_clonal_hematopoiesis/HC_vcf_link/','DB8_48s_HC_SNVs_merged.txt.inc8_Blood_3s_merged.less5.readinfo','8_Blood_3s_merged',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/10_clonal_hematopoiesis/HC_vcf_link/','DB9_36s_HC_SNV_merged.txt.inc9_Blood_3s_merged.less5.readinfo','9_Blood_3s_merged',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/10_clonal_hematopoiesis/HC_vcf_link/','DB10_40s_HC_SNV_merged.txt.inc10_Blood_3s_merged.less5.readinfo','10_Blood_3s_merged'
)


#Filter1 with readinfo
if(F){
  path_tbl2 <- path_tbl2 %>% mutate(path1 = paste0(folder_path, file_name1)) #SNV
  path_tbl2 <- path_tbl2 %>% mutate(path1 = gsub("SNV", "indel",paste0(folder_path, file_name1))) #indel
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    print(this_db)
    dt <- read_tsv(path_tbl2$path1[path_tbl2$deadbody == this_db], col_types = cols(`#CHROM`='c'))
    m_dt <- dt %>% separate(`ref_minMQ;med;max`, c("ref_minMQ","ref_medMQ","ref_maxMQ"), convert= T, sep =';') %>%
      separate(`var_minMQ;med;max`, c("var_minMQ","var_medMQ","var_maxMQ"), covert = T, sep =';') %>%
      separate(`ref_med;meanBQ`, c("ref_medBQ","ref_meanBQ"), convert = T, sep=';') %>%
      separate(`var_med;meanBQ`, c("var_medBQ","var_meanBQ"), convert =T, sep=';') %>%
      separate(`var_minMismatch;med;max`, c("var_minMismatch","var_medMismatch","var_maxMisMatch"), convert = T, sep=';') %>%
      separate(`refClip;pct;varClip;pct`, c("refClipN","refClipPct","varClipN","varClipPct"), convert = T, sep=';') %>%
      separate(`var_LocaLeftMin;Med;Max`, c("var_LocaLeftMin","var_LocaLeftMed","var_LocaLeftMax"), convert = T, sep=';') %>%
      separate(`var_LocaRightMin;Med;Max`, c("var_LocaRightMin","var_LocaRightMed","var_LocaRightMax"), convert = T, sep=';') %>%
      mutate_at(vars(starts_with("var_")), as.numeric) %>%
      mutate(VAF = (var_readN)/(var_readN + ref_readN)) %>%
      mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
    dim(m_dt)
    subdt <- m_dt %>% filter(var_medMQ >= 25 & 
                               ref_medMQ - var_medMQ < 10 & 
                               var_meanBQ >= 20 & 
                               #NP_maxVAF < 1 & 
                               var_readN > 2 &
                               ref_minMQ > 0 & 
                               varClipPct < 95 &
                               #!(var_id %in% false_list[[this_db]]) & 
                               var_medMismatch < 5)
    print(nrow(subdt))
    write_tsv(subdt,paste0(path_tbl2$path1[path_tbl2$deadbody == this_db],'.fi'))
  }
}


#Filter2 with other sample info
if(F){
  path_tbl2 <- path_tbl2 %>% mutate(path2 = paste0(folder_path, file_name1,'.fi.readc'))
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    print(this_db)
    dt <- read_tsv(path_tbl2$path2[path_tbl2$deadbody == this_db], col_types = cols(`#CHROM`='c'))
    colnames(dt)[(ncol(dt)-9):ncol(dt)] <- paste0('S',c(1:10))
    m_dt <- dt %>% 
      separate(S1, c('S1_ref','S1_var','S1_vaf'),convert = T, sep=';') %>%
      separate(S2, c('S2_ref','S2_var','S2_vaf'),convert = T, sep=';') %>%
      separate(S3, c('S3_ref','S3_var','S3_vaf'),convert = T, sep=';') %>%
      separate(S4, c('S4_ref','S4_var','S4_vaf'),convert = T, sep=';') %>%
      separate(S5, c('S5_ref','S5_var','S5_vaf'),convert = T, sep=';') %>%
      separate(S6, c('S6_ref','S6_var','S6_vaf'),convert = T, sep=';') %>%
      separate(S7, c('S7_ref','S7_var','S7_vaf'),convert = T, sep=';') %>%
      separate(S8, c('S8_ref','S8_var','S8_vaf'),convert = T, sep=';') %>%
      separate(S9, c('S9_ref','S9_var','S9_vaf'),convert = T, sep=';') %>%
      separate(S10, c('S10_ref','S10_var','S10_vaf'),convert = T, sep=';') %>%
      mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
    
    m_dt <- left_join(m_dt, m_dt %>% select(var_id, ends_with('_vaf')) %>% mutate(n_posi = apply(.[2:11], 1, function(x) sum(x>0)),
                                                                                  n_less25 = apply(.[2:11], 1, function(x) sum(x<25 & x >0))) %>% select(var_id, n_posi, n_less25))
    
    sub_dt <- m_dt %>% filter(n_posi < 3)
    print(nrow(sub_dt))
    write_tsv(sub_dt, paste0(path_tbl2$path2[path_tbl2$deadbody == this_db], '.fi'))
  }
}

#heatmap generation for final manual filter
if(F){
  path_tbl2 <- path_tbl2 %>% mutate(n_sample = apply(.["file_name1"], 1, function(x) unlist(strsplit(x, '_'))[2])) %>%
    mutate(path3 = paste0(path2, '.fi.', n_sample,'Call'))
  path_tbl2 %>% View()
  library(ComplexHeatmap)
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    print(this_db)
    dt <- read_tsv(path_tbl2$path3[path_tbl$deadbody == this_db],
                   col_types = cols(`#CHROM`='c'))
    mx <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) %>%
      select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
    pdf(paste0(path_tbl2$path3[path_tbl$deadbody == this_db],'.heatmap.pdf'))
    dt <- left_join(dt, dt %>% select(var_id,ends_with('vafpct')) %>% mutate(n_posi_total = apply(.[2:ncol(.)], 1, function(x) sum(x>0))) %>%
                      select(var_id, n_posi_total))
    subdt <- dt %>% filter(n_posi_total ==1)
    print(nrow(subdt))
    print(Heatmap(mx))
    dev.off()
  }
}

#final list generation
if(F){
  path_tbl2 <- path_tbl2 %>% mutate(path4 = paste0(path2, '.fi'))
  pass_list <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/clonal_hematopoiesis/CH_passed_variants_list.txt') %>% pull(var_id)
  merged_tbl <- tibble()
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    print(this_db)
    dt <- read_tsv(path_tbl2$path4[path_tbl$deadbody == this_db],
                   col_types = cols(`#CHROM`='c'))
    
    subdt <- dt %>% filter(var_id %in% pass_list)
    merged_tbl <- bind_rows(merged_tbl, subdt %>% mutate(deadbody = this_db))
    print(nrow(subdt))
    write_tsv(subdt,paste0(path_tbl2$path4[path_tbl$deadbody == this_db],'.fi'))
  }
  write_tsv(merged_tbl,'/home/users/sypark/00_Project/06_LineageTracing/merged_data/clonal_hematopoiesis/CH_final_list.txt')
}




#somatic call filter
nc_change <- c("TG","TC","TA","CA","CG","CT","CT","CG","CA","TA","TC","TG")
names(nc_change) <- c("AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG")



snv_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/01_BAM/6_blood_LINFb1_strelka2/results/variants/somatic.snvs.vcf.pass.readinfo.readc.rasm',
               comment='##', col_types =cols(`#CHROM`='c'))
snv_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db9/10_clonal_hematopoiesis/02_strelka2/9_blood_ARUL1400F8_strelka2/results/variants/somatic.snvs.vcf.pass.readinfo.readc.rasm',
                   comment='##', col_types =cols(`#CHROM`='c'))

colnames(snv_dt)
snv_dt <- snv_dt %>% mutate(nc_type = nc_change[paste0(REF,ALT)]) %>%
  mutate(vaf = var_readN/(ref_readN + var_readN)) %>%
  separate(`ref_minMQ;med;max`, c('ref_minMQ','ref_medMQ','ref_maxMQ'), sep=';', convert=T) %>%
  separate(`var_minMQ;med;max`, c('var_minMQ','var_medMQ','var_maxMQ'), sep=';', convert=T) %>%
  separate(`ref_med;meanBQ`, c('ref_medBQ','ref_meanBQ'), sep=';', convert=T) %>%
  separate(`var_med;meanBQ`, c('var_medBQ','var_meanBQ'), sep=';', convert=T) %>%
  separate(`var_LocaLeftMin;Med;Max`,c('var_LocaLeftMin','var_LocaLeftMed','var_LocaLeftMax'), sep=';', convert=T) %>%
  separate(`var_LocaRightMin;Med;Max`,c('var_LocaRightMin','var_LocaRightMed','var_LocaRightMax'), sep=';', convert=T) %>%
  separate(`refClip;pct;varClip;pct`, c('refClip','refClip_pct','varClip','varClip_pct'), sep=';', convert=T) %>%
  separate(`refIns;pct;varIns;pct`, c('refIns','refIns_pct','varIns','varIns_pct'), sep=';', convert=T) %>%
  separate(`refDel;pct;varDel;pct`, c('refDel','refDel_pct','varDel','varDel_pct'), sep=';', convert=T) %>%
  separate(`ref_minMismatch;med;max`, c('ref_minMismatch','ref_medMismatch','ref_maxMismatch'), sep=';', convert=T) %>%
  separate(`var_minMismatch;med;max`, c('var_minMismatch','var_medMismatch','var_maxMismatch'), sep=';', convert=T) %>%
  separate(`pairN_ref_readN;var_readN;vaf%`, c('pairN_ref_readN','pairN_var_readN','pairN_vaf_pct'), sep=';', convert=T) %>%
  separate("repeat_unit;ref_repeat_count;lt_seq_size;rt_seq_size;t_ref;t_var;t_unknown;n_other;n_consen;t_ref_cor;t_var_cor;t_ukn_cor", c('repeat_unit','ref_repeat_count','lt_seq_size','rt_seq_size','t_ref','t_var','t_ukn','n_other','n_consen','t_ref_cor','t_var_cor','t_ukn_cor'), sep=';', convert=T)

nrow(snv_dt)
chr_list <- c(1:22,'X','Y')
f_snv_dt <- snv_dt %>% filter(`#CHROM` %in% chr_list & ref_medMQ  - var_medMQ < 10 & pairN_vaf_pct < 1 & n_consen ==0) 
nrow(f_snv_dt)
#f_snv_dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/strelka2_pass_RI_fi.tsv')
f_snv_dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db9/10_clonal_hematopoiesis/02_strelka2/9_blood_ARUL1400F8_strelka2/results/variants/somatic.snvs.vcf.pass.readinfo.readc.rasm.ri_fi')
ggplot(f_snv_dt, aes(x=vaf))+
  geom_histogram()
nrow(f_snv_dt)
ggplot(f_snv_dt, aes(x=var_readN, y=vaf))+
  geom_jitter(aes(color=nc_type), alpha=0.5)+
  geom_abline(slope=-0.01, intercept=0.15)



#after readinfo filter and mpu
mpu_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/02_strelka2/6_blood_LINFb1_strelka2_pass_RI_fi.tsv.10sCall',
                   comment='##', col_types =cols(`#CHROM`='c'))
mpu_dt <- mpu_dt %>% mutate(vafpct = 100*vaf, 
                  var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) 
colnames(mpu_dt)
s_mpu_dt <- mpu_dt %>%select(var_id, vafpct, ends_with('vafpct'))
mx <- s_mpu_dt %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
library(pheatmap)
pheatmap(mx)

s_mpu_dt$n_m0 <- rowSums(mx > 0)
f_mpu_dt <- s_mpu_dt %>% filter(n_m0 <= 2) %>% select(-n_m0)
f_mx <- f_mpu_dt %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
pheatmap(f_mx)
nrow(f_mx)
f_ids <- rownames(f_mx)

f_mpu_dt <- mpu_dt %>% filter(var_id %in% f_ids)
f_mpu_dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/02_strelka2/6_blood_LINFb1_strelka2_pass_RI_fi.tsv.10sCall.fi')
ggplot(f_mpu_dt, aes(x=vafpct))+
  geom_histogram()

ggplot(f_mpu_dt, aes(x=var_readN, y=vaf))+
  geom_jitter(aes(color=nc_type), alpha=0.5)

ggplot(f_mpu_dt)+
  geom_histogram(aes(x=vaf, fill=nc_type), alpha=0.5)

f_mpu_dt %>% filter(vaf>0.14 & vaf <0.2)  

#after additional pileup
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/02_strelka2/6_blood_LINFb1_strelka2_pass_RI_fi.tsv.10sCall.fi.50sCall',
               col_types= cols(`#CHROM`='c'))
s_dt <- dt %>%select(var_id, vafpct, ends_with('vafpct'))
s_mx <- s_dt %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
library(pheatmap)
pheatmap(s_mx)
dt$n_m0 <- rowSums(s_mx > 0)
f_dt <- dt %>% filter(n_m0 <= 2) %>% select(-n_m0)
f_mx <- dt %>%select(var_id, vafpct, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
pheatmap(f_mx)
nrow(f_mx)
#f_dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/02_strelka2/6_blood_LINFb1_strelka2_pass_RI_fi.tsv.10sCall.fi.50sCall.fi')
ggplot(f_dt, aes(x=var_readN, y=vaf))+
  geom_jitter(aes(color=nc_type), alpha=0.5)+
  geom_abline(slope=-0.01, intercept=0.15)

ggplot(f_dt, aes(x=vaf))+
  geom_histogram()

f_dt %>% filter(vaf > var_readN*(-0.01)+0.15) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/13_clonal_hematopoiesis/02_strelka2/tmp_fi.txt')
