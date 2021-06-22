
###########################<Common processing>##############################################
#load library
library(tidyverse)
library(ggsci)
library(ggtree)
library(ggrepel)
library(cowplot)
library(scales)
library(grid)
library(gridExtra)
library(png)

#load dataset
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_200918.txt')

# path assign
nwk_tbl <- tribble(
  ~deadbody, ~nwk_path,
  "DB2", '/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/DB2_Lineage_count_table.txt.nwk', 
  "DB3", '/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/DB3_Lineage_count_table.txt.nwk',
  "DB5",'/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/DB5_Lineage_count_table.txt.nwk',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/DB6_Lineage_count_table.txt.nwk',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/DB8_Lineage_count_table.txt.nwk',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/DB9_Lineage_count_table.txt.nwk',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/DB10_Lineage_count_table.txt.nwk'
)
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

#edit meta_dt
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit lng_dt
lng_dt$deadbody <- factor(lng_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit DB_dt
DB_dt$deadbody <- factor(DB_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#color setting

#show_col(pal_npg('nrc')(10))
db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)
cell_pal <- pal_npg('nrc')(6)
names(cell_pal) <- unique(meta_dt$Cell_type)
source2_pal <- pal_d3('category10')(6)
names(source2_pal) <- c('HN','lt_LE','lt_UE','rt_LE','trunk','rt_UE')
timing_pal <- pal_npg('nrc')(10)[c(1,3,4)]
names(timing_pal) <- c('early','early_to_late','late')
pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
sig_pal <- pal_combi[c(9,11,1,12,5,17,13)]
names(sig_pal) <- c('v3_1','v3_5','v3_7a','v3_7b','v3_7c','v3_7d','v3_18')


#Fx source
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/merged_analysis_functions.R')

#theme setting
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), strip.text = element_text(size=15), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))






##############################################<<CNV>>>######################################
#CNV_analysis
#making list of CNV more than 10Mbps
if(F){
  file_list = list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_seg/', pattern= 'segments', full.names = T)
  cyto_dt <- read_tsv("/home/users/data/02_annotation/02_annovar/humandb_hg19/hg19_cytoBand.txt", col_names = c('chrom', 'start.pos','end.pos','cyto_id','cyto_stain'))
  acen_dt <- cyto_dt %>% filter(cyto_stain == 'acen' | (cyto_stain == 'gvar' & chrom != 'chrY')) %>% group_by(chrom) %>% summarise(cen_start = min(start.pos), cen_end = max(end.pos))
  acen_dt$chrom <- gsub('chr','',acen_dt$chrom)
  acen_dt$cen_end[acen_dt$chrom == '1'] <- 150000000  #manual correction
  cyto_dt$chrom <- gsub('chr','',cyto_dt$chrom)
  f_cyto_dt <- left_join(cyto_dt %>% filter(!cyto_stain %in% c('acen','gvar','stalk')), acen_dt) %>% mutate(arm = ifelse(end.pos <= cen_start, 'p', 'q'))
  f_cyto_dt <- f_cyto_dt %>% group_by(chrom, arm) %>% summarise(start.pos = min(start.pos), end.pos = max(end.pos)) %>% mutate(dist = end.pos - start.pos)
  possible_ps <- f_cyto_dt %>% filter(arm == 'p') %>% pull(chrom) %>% unique()
  merged_tbl <- tibble()
  for (i in 1:length(file_list)){
    print(i)
    dt <- read_tsv(file_list[i], col_types = cols(`#CHROM`='c'))
    t_db <- paste0('DB',unlist(strsplit(unlist(strsplit(file_list[i],'//'))[2],'_'))[1])
    t_id <- unlist(strsplit(unlist(strsplit(file_list[i],'//'))[2],'\\.'))[1]
    print(t_id)
    t_gender <- unique(meta_dt$gender[meta_dt$deadbody == t_db])
    dt <- dt %>% rename(chrom= `#CHROM`)
    total_dt <- dt %>% group_by(chrom) %>% summarise(chrom_start = min(start_pos), chrom_end = max(end_pos)) %>% mutate(total_dist = chrom_end - chrom_start)
    if( t_gender == 'F'){
      amp_dt <- dt %>% filter(majCN > 1) %>% mutate(dist = end_pos - start_pos) %>% mutate(type = 'amp')
      del_dt <- dt %>% filter(minCN < 1) %>% mutate(dist = end_pos - start_pos) %>% mutate(type = 'del')
      final_dt <- left_join(bind_rows(amp_dt, del_dt), total_dt) %>% mutate(ratio = dist/total_dist, deadbody = t_db, sample_id = t_id)
    } else if (t_gender == 'M'){ 
      amp_dt <- dt %>% filter(majCN > 1) %>% mutate(dist = end_pos - start_pos) %>% mutate(type = 'amp')
      del_dt <- dt %>% filter((chrom %in% 1:22 & minCN < 1) | chrom %in% c('X','Y') & majCN <1) %>% mutate(dist = end_pos - start_pos) %>% mutate(type = 'del')
      final_dt <- left_join(bind_rows(amp_dt, del_dt), total_dt) %>% mutate(ratio = dist/total_dist, deadbody = t_db, sample_id = t_id)
    }
    merged_tbl <- bind_rows(merged_tbl, final_dt)
  }
  merged_tbl <- left_join(merged_tbl,acen_dt)
  residual = 3000000
  near_co = 5000000
  
  f_merged_tbl <- merged_tbl %>% filter(dist > 10000000)  #filter dist > 10Mbp
  m_merged_tbl <- f_merged_tbl %>% mutate(group1 = ifelse(ratio > 0.95, 'whole_chrom', 
                                                          ifelse(start_pos >= (cen_start -residual) & start_pos <= (cen_end + residual) & abs(end_pos - chrom_end) <= residual, 'q_arm',
                                                                 ifelse(abs(start_pos - chrom_start) <= residual & end_pos >= (cen_start - residual) & end_pos <= (cen_end + residual), 'p_arm', 
                                                                        ifelse(start_pos >= cen_start & end_pos <= cen_end, 'intra_cent',
                                                                               ifelse(abs(start_pos - cen_start) < near_co | abs(end_pos - cen_end) < near_co, 'near_cent',
                                                                                      'intra_chrom'
                                                                               ))))))
  m_merged_tbl <- m_merged_tbl %>% filter(chrom %in% possible_ps | group1 != 'p_arm')  #filter false segments
  m_merged_tbl <- m_merged_tbl %>% filter(group1 != 'near_cent') #filter near centromeric false segments
  m_merged_tbl <- m_merged_tbl %>% mutate(cnv_id = ifelse(group1 %in% c('whole_chrom','q_arm'), paste(group1, type, chrom, sep='_'), paste(group1,type,chrom,start_pos,end_pos,sep='_') ))
  saveRDS(m_merged_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.rds')
}

if(F){
  m_merged_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.rds')
  
  #Manual subclonality check
  m_merged_tbl %>% select(-chrom_start, -chrom_end, -total_dist, -ratio) %>% View()
  if(F){
    subclonal_dt %>% saveRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.subclonalCN.rds')
  }
  
  #remove some subclonal CNVs
  subclonal_dt <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.subclonalCN.rds')
  subclonal_dt <- subclonal_dt %>% mutate(subcl = 'subclonal')
  subclonal_dt <- subclonal_dt[-4,]
  m_merged_tbl <- left_join(m_merged_tbl, subclonal_dt %>% select(chrom, start_pos, end_pos, sample_id, subcl))
  m_merged_tbl <- m_merged_tbl %>% filter(is.na(subcl)) %>% select(-subcl)
  if(F){
    saveRDS(m_merged_tbl,'/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.rmSC.rds')
  }
  
}


m_merged_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.rmSC.rds')
#list of DB3 CNVs (Fig 4i)
m_merged_tbl %>% filter(deadbody == 'DB3') %>% group_by(chrom) %>% count()

#CNV analysis
cnv_tbl <- m_merged_tbl %>% select(sample_id, cnv_id) %>% mutate(presence='Y') %>% spread(key=cnv_id, value=presence)
cnv_tbl[is.na(cnv_tbl)] <- 'N'
m_lng_dt <- #join lng_dt and cnv_tbl
  left_join(lng_dt,left_join(cnv_tbl, meta_dt %>% mutate(lineage_id2 = paste(deadbody, lineage_id_cindel, sep="_")) %>% select(lineage_id2, sample_id)),by='lineage_id2')
node_list <- m_lng_dt %>% filter(is.na(sample_id) == T) %>% .$lineage_id2
tail(node_list)
for (node_id in node_list){   #write CNV info per lineage_id2, common CNV in descendants -> CNV of ancestor
  print(node_id)
  for(col_id in setdiff(colnames(cnv_tbl), c('sample_id'))){
    ex_values <- m_lng_dt %>% filter(grepl(paste0(node_id,'-'), lineage_id2)==T & is.na(get(col_id)) == F) %>% .[[col_id]] %>% unique()
    print(ex_values)
    if (length(ex_values)==1){
      m_lng_dt[which(m_lng_dt$lineage_id2 == node_id),][[col_id]] <- ex_values
    } else{
      m_lng_dt[which(m_lng_dt$lineage_id2 == node_id),][[col_id]] <- 'N'
    }
  }
}
#count CNV event in early_to_late lineages
fm_lng_dt <- m_lng_dt %>% filter(time_of_lineage == 'early_to_late')
tmp_tbl <- fm_lng_dt %>% select(setdiff(colnames(cnv_tbl), c('sample_id'))) 
mx <- tmp_tbl %>% as.matrix()
fm_lng_dt$anyCNV <- rowSums(mx == 'Y')
mx <- tmp_tbl %>% select(starts_with('whole_chrom_amp')) %>% as.matrix()
fm_lng_dt$anyWCA <- rowSums(mx == 'Y')
mx <- tmp_tbl %>% select(starts_with('whole_chrom_del')) %>% as.matrix()
fm_lng_dt$anyWCD <- rowSums(mx == 'Y')
mx <- tmp_tbl %>% select(starts_with('q_arm_amp')) %>% as.matrix()
fm_lng_dt$anyQAA <- rowSums(mx == 'Y')
mx <- tmp_tbl %>% select(starts_with('q_arm_del')) %>% as.matrix()
fm_lng_dt$anyQAD <- rowSums(mx == 'Y')
mx <- tmp_tbl %>% select(starts_with('intra_chrom_amp')) %>% as.matrix()
fm_lng_dt$anyICA <- rowSums(mx == 'Y')
mx <- tmp_tbl %>% select(starts_with('intra_chrom_del')) %>% as.matrix()
fm_lng_dt$anyICD <- rowSums(mx == 'Y')

#
tmp1 <- fm_lng_dt %>% select(deadbody, lineage_id, samples, starts_with('any'))
tmp1$sample_n <- map_int(tmp1$samples, function(x) length(unlist(strsplit(x, ';'))))
tmp1

#No. of CNV events per deadbody
type_tbl <- tibble(deadbody = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
ct_tbl <- fm_lng_dt %>% group_by(deadbody, anyCNV) %>% count() %>% mutate(anyCNVpre = ifelse(anyCNV == 0,'N','Y')) %>% 
  group_by(deadbody, anyCNVpre) %>% summarise(anyCNVpre_n= sum(n))
total_cts <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(total_ct = n())
ct_tbl <- left_join(ct_tbl,total_cts) %>% filter(anyCNVpre == 'Y') %>% mutate(anyCNVprop = anyCNVpre_n/total_ct) %>% select(deadbody, total_ct, anyCNVpre_n, anyCNVprop)
type_tbl <- left_join(type_tbl, ct_tbl)
ct_tbl <- fm_lng_dt %>% group_by(deadbody, anyWCA) %>% count() %>% mutate(anyWCApre = ifelse(anyWCA == 0,'N','Y')) %>% 
  group_by(deadbody, anyWCApre) %>% summarise(anyWCApre_n= sum(n))
total_cts <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(total_ct = n())
ct_tbl <- left_join(ct_tbl,total_cts) %>% filter(anyWCApre == 'Y') %>% mutate(anyWCAprop = anyWCApre_n/total_ct) %>% select(deadbody, anyWCApre_n, anyWCAprop)
type_tbl <- left_join(type_tbl, ct_tbl)
ct_tbl <- fm_lng_dt %>% group_by(deadbody, anyWCD) %>% count() %>% mutate(anyWCDpre = ifelse(anyWCD == 0,'N','Y')) %>% 
  group_by(deadbody, anyWCDpre) %>% summarise(anyWCDpre_n= sum(n))
total_cts <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(total_ct = n())
ct_tbl <- left_join(ct_tbl,total_cts) %>% filter(anyWCDpre == 'Y') %>% mutate(anyWCDprop = anyWCDpre_n/total_ct) %>% select(deadbody, anyWCDpre_n, anyWCDprop)
type_tbl <- left_join(type_tbl, ct_tbl)
ct_tbl <- fm_lng_dt %>% group_by(deadbody, anyQAA) %>% count() %>% mutate(anyQAApre = ifelse(anyQAA == 0,'N','Y')) %>% 
  group_by(deadbody, anyQAApre) %>% summarise(anyQAApre_n= sum(n))
total_cts <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(total_ct = n())
ct_tbl <- left_join(ct_tbl,total_cts) %>% filter(anyQAApre == 'Y') %>% mutate(anyQAAprop = anyQAApre_n/total_ct) %>% select(deadbody, anyQAApre_n, anyQAAprop)
type_tbl <- left_join(type_tbl, ct_tbl)
ct_tbl <- fm_lng_dt %>% group_by(deadbody, anyQAD) %>% count() %>% mutate(anyQADpre = ifelse(anyQAD == 0,'N','Y')) %>% 
  group_by(deadbody, anyQADpre) %>% summarise(anyQADpre_n= sum(n))
total_cts <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(total_ct = n())
ct_tbl <- left_join(ct_tbl,total_cts) %>% filter(anyQADpre == 'Y') %>% mutate(anyQADprop = anyQADpre_n/total_ct) %>% select(deadbody, anyQADpre_n, anyQADprop)
type_tbl <- left_join(type_tbl, ct_tbl)
ct_tbl <- fm_lng_dt %>% group_by(deadbody, anyICA) %>% count() %>% mutate(anyICApre = ifelse(anyICA == 0,'N','Y')) %>% 
  group_by(deadbody, anyICApre) %>% summarise(anyICApre_n= sum(n))
total_cts <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(total_ct = n())
ct_tbl <- left_join(ct_tbl,total_cts) %>% filter(anyICApre == 'Y') %>% mutate(anyICAprop = anyICApre_n/total_ct) %>% select(deadbody, anyICApre_n, anyICAprop)
type_tbl <- left_join(type_tbl, ct_tbl)
ct_tbl <- fm_lng_dt %>% group_by(deadbody, anyICD) %>% count() %>% mutate(anyICDpre = ifelse(anyICD == 0,'N','Y')) %>% 
  group_by(deadbody, anyICDpre) %>% summarise(anyICDpre_n= sum(n))
total_cts <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(total_ct = n())
ct_tbl <- left_join(ct_tbl,total_cts) %>% filter(anyICDpre == 'Y') %>% mutate(anyICDprop = anyICDpre_n/total_ct) %>% select(deadbody, anyICDpre_n, anyICDprop)
type_tbl <- left_join(type_tbl, ct_tbl)
type_tbl[is.na(type_tbl)] <-0

#Proportion of samples with CNV per deadbody
type_tbl$anyCNVprop %>% summary()
type_tbl$deadbody <- factor(type_tbl$deadbody, levels = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
ggplot(type_tbl, aes(x=deadbody, y=anyCNVprop, fill= deadbody))+
  geom_bar(stat="identity") +
  scale_fill_manual(values= db_pal) +
  ylab("proportion") +
  ggtitle("Proportion of samles with any CNV (>10Mbp)")+
  theme_syp+theme(axis.title.x = element_blank())

#Proportion of samples with CNV per CNV type per deadbody
m_type_tbl <- type_tbl %>% select(-anyCNVprop) %>% gather(-deadbody, key='type', value='proportion')
ggplot(m_type_tbl, aes(x=type, y=proportion, fill=type))+
  geom_bar(stat="identity", position = "dodge")+
  facet_wrap(~deadbody)+
  ylab("Proportion of samples with CNV")+ xlab("CNV type")+
  theme(axis.text.x = element_text(angle=45, hjust=1))

#event number per CNVtype (EDFig 5b)
event_tbl <- fm_lng_dt %>% select(deadbody, lineage_id2, starts_with('any')) %>% group_by(deadbody) %>% summarise(CNVevent_n = sum(anyCNV), WCAevent_n = sum(anyWCA), WCDevent_n = sum(anyWCD),
                                                                                                                  QAAevent_n = sum(anyQAA), QADevent_n = sum(anyQAD), ICAevent_n = sum(anyICA),
                                                                                                                  ICDevent_n = sum(anyICD))
m_event_tbl <- event_tbl %>% select(-CNVevent_n) %>% gather(-deadbody, key="type", value="count")
type2_v <- c('WCA','WCD','AA','AD','ICA','ICD')
names(type2_v) <- m_event_tbl$type %>% unique()
m_event_tbl$type2 <- type2_v[m_event_tbl$type]
m_event_tbl$type2 <- factor(m_event_tbl$type2, levels = c('WCA','WCD','AA','AD','ICA','ICD'))
ggplot(m_event_tbl, aes(x=type2, y=count))+
  geom_bar(aes(fill = deadbody), stat="identity")+
  scale_fill_manual(values = db_pal)+
  xlab("CNV type")+ylab("Count")+
  theme_syp+theme(axis.text.x = element_text(angle=45, hjust=1))
m_event_tbl$count %>% sum()
m_event_tbl$type %>% unique()
m_event_tbl %>% filter(type %in% c('WCAevent_n','WCDevent_n')) %>% pull(count)  %>% sum()

#event number per chromosome (EDFig 5c)
colnames(fm_lng_dt)
tmp_tbl <- fm_lng_dt %>% select(deadbody, lineage_id2, starts_with('whole_chrom_'), starts_with('q_arm_'), starts_with('intra_chrom_')) 
tmp_tbl[tmp_tbl == 'N'] <- 0
tmp_tbl[tmp_tbl == 'Y'] <- 1
tmp_tbl <- tmp_tbl %>% mutate_at(vars(-deadbody, -lineage_id2), as.numeric)
tmp_tbl <- tmp_tbl %>% group_by(deadbody) %>% summarise_at(vars( -lineage_id2), list(sum=sum))
m_tmp_tbl <- tmp_tbl %>% gather(-deadbody, key=type, value=count)
m_tmp_tbl <- m_tmp_tbl %>% mutate(chrom = apply(.["type"], 1, function(x) unlist(strsplit(x, '_'))[4]))  %>%
  mutate(cnv_type = apply(.["type"], 1, function(x) paste(unlist(strsplit(x, '_'))[1:3], collapse="_")))
m_tmp_tbl <- bind_rows(m_tmp_tbl, tibble(deadbody =c('DB2','DB2'),type =c('.','.'),count = c(0,0),chrom = c("9","19"),cnv_type = c('intra_chrom_amp', 'intra_chrom_amp'))) #add blank row for chrom9,19
m_tmp_tbl$chrom = factor(m_tmp_tbl$chrom, levels =c(1:22,'X','Y'))
m_tmp_tbl$cnv_type %>% unique()
fm_tmp_tbl <- m_tmp_tbl %>% filter(cnv_type %in% c("whole_chrom_amp","whole_chrom_del"))
fm_tmp_tbl$count %>% sum()
fm_tmp_tbl$count[fm_tmp_tbl$chrom %in% c('X','Y')] %>% sum()
remainder = setdiff(c(1:22,'X','Y'),fm_tmp_tbl$chrom %>% unique())
fm_tmp_tbl <- bind_rows(fm_tmp_tbl, tibble(deadbody = 'DB2', type = 'none',count = 0, chrom = remainder, cnv_type = 'none'))
fm_tmp_tbl$chrom = factor(fm_tmp_tbl$chrom, levels =c(1:22,'X','Y'))
ggplot(fm_tmp_tbl, aes(x=chrom, y=count, fill=deadbody))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = db_pal)+
  xlab("Chromosome")+ylab("No. of events")+
  scale_y_continuous(limits = c(0,12), breaks = seq(0,12,2))+
  theme_syp


#proportion of samples with sex chromosomal losses
tmp_tbl2 <- left_join(m_tmp_tbl%>% filter(chrom %in% c('X','Y') & cnv_type == 'whole_chrom_del') %>% group_by(deadbody) %>% summarise(count=sum(count)), total_cts)
tmp_tbl2 <- tmp_tbl2 %>% mutate(ratio = count/total_ct)
tmp_tbl2$deadbody = factor(tmp_tbl2$deadbody, levels = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
ggplot(tmp_tbl2, aes(x=deadbody, y=ratio, fill= deadbody))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = db_pal)+
  theme_syp+
  ggtitle('Proportion of samples with sex chromosome loss')


m_merged_tbl



#####################CNV timing###############
#plot mutCN per AMP
if(F){
  m_merged_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.rmSC.rds')
  f_merged_tbl <- m_merged_tbl %>% filter(type == 'amp' & dist > 50000000)
  early_v = c()
  l_early_v = c()
  h_early_v = c()
  late_v = c()
  l_late_v = c()
  h_late_v = c()
  adj_ratio_v = c()
  l_adj_ratio_v = c()
  h_adj_ratio_v = c()
  for (i in 1:nrow(f_merged_tbl)){
    print(i)
    t_id <- f_merged_tbl[i,"sample_id"] %>% as.character()
    t_db <- f_merged_tbl[i,"deadbody"] %>% as.character()
    t_path1 <- paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/timing_analysis/amp_timing/',t_id,'.edited.vcf.anv.seqzcn.bino1-4P')
    t_chrom <- f_merged_tbl[i,"chrom"] %>% as.character()
    t_start <- f_merged_tbl[i,"start_pos"] %>% as.numeric()
    t_end <- f_merged_tbl[i,"end_pos"] %>% as.numeric()
    t_majCN <- f_merged_tbl[i,"majCN"] %>% as.numeric()
    t_minCN <- f_merged_tbl[i,"minCN"] %>% as.numeric()
    t_dist <- paste0(round((t_end - t_start)/1000000,1),'Mbps')
    
    dt <- read_tsv(t_path1, comment= '##', col_types = cols(`#CHROM`='c', pEarly ='n',`pLate/minor` ='n')) %>% filter(POS >= t_start & POS <= t_end & `#CHROM` == t_chrom)
    dt <- dt %>% separate(get(t_id), c('RD','AD','VAF'), convert = T, sep=':') %>% mutate(sample_id = t_id, TCF = 1, mploidy=2) %>% 
      rename(pLate_minor = `pLate/minor`)
    dt <- dt %>% mutate(mutCN = VAF*0.01*totCN*(totCN*TCF + (1-TCF)*2)/totCN/TCF, earlyV = pEarly*(1-pEarly), lateV = pLate_minor * (1-pLate_minor))
    early_n <- sum(dt$pEarly)
    l_early_n <- early_n - 1.96*sqrt(sum(dt$earlyV))
    h_early_n <- early_n + 1.96*sqrt(sum(dt$earlyV))
    late_n <- sum(dt$pLate_minor)
    l_late_n <- late_n - 1.96*sqrt(sum(dt$lateV))
    h_late_n <- late_n + 1.96*sqrt(sum(dt$lateV))
    if(t_majCN > t_minCN & t_minCN > 0){
      adj_early_n = early_n*2
      adj_l_early_n = l_early_n *2
      adj_h_early_n = h_early_n*2
      adj_late_n = max(0,late_n - early_n)
      adj_l_late_n = max(0, l_late_n - h_early_n)
      adj_h_late_n = max(0, h_late_n - l_early_n)
    } else{  
      adj_early_n = early_n
      adj_l_early_n = l_early_n
      adj_h_early_n = h_early_n
      adj_late_n = late_n
      adj_l_late_n = l_late_n
      adj_h_late_n = h_late_n
    } 
    adj_ratio = round(adj_early_n/(adj_early_n + adj_late_n),2)
    l_adj_ratio = round(adj_l_early_n/(adj_l_early_n + adj_h_late_n),2)
    h_adj_ratio = round(adj_h_early_n/(adj_h_early_n + adj_l_late_n),2)
    
    early_v = c(early_v, early_n)
    l_early_v = c(l_early_v, l_early_n)
    h_early_v = c(h_early_v, h_early_n)
    late_v = c(late_v, late_n)
    l_late_v = c(l_late_v, l_late_n)
    h_late_v = c(h_late_v, h_late_n)
    adj_ratio_v = c(adj_ratio_v, adj_ratio)
    l_adj_ratio_v = c(l_adj_ratio_v, l_adj_ratio)
    h_adj_ratio_v = c(h_adj_ratio_v, h_adj_ratio)
  }
  
  #add info acquired during loop
  f_merged_tbl$early_n <- early_v
  f_merged_tbl$early_n_low <- l_early_v
  f_merged_tbl$early_n_hi <- h_early_v
  f_merged_tbl$late_minor_n <- late_v
  f_merged_tbl$late_minor_n_low <- l_late_v
  f_merged_tbl$late_minor_n_hi <- h_late_v
  f_merged_tbl$adj_ratio <- adj_ratio_v
  f_merged_tbl$adj_ratio_low <- l_adj_ratio_v
  f_merged_tbl$adj_ratio_hi <- h_adj_ratio_v
  
  #manual averaging of late diverging samples
  mf_merged_tbl <- left_join(f_merged_tbl,meta_dt %>% select(deadbody, age) %>% unique()) %>% 
    mutate(CNV_time = adj_ratio *age, CNV_time_low = adj_ratio_low * age, CNV_time_hi = adj_ratio_hi * age) 
  mf_merged_tbl <- mf_merged_tbl %>% mutate(sample_id2 = ifelse(sample_id %in% c('3_LP_Fb_001G9','3_LP_Fb_001B12a','3_LP_Fb_001B3a'), 'DB3_L1-1',sample_id)) 
  mf_merged_tbl <- mf_merged_tbl %>% group_by(deadbody, age, group1,cnv_id,sample_id2) %>% 
    summarise(CNV_time = mean(CNV_time), CNV_time_low = mean(CNV_time_low), CNV_time_hi = mean(CNV_time_hi),
              adj_ratio = mean(adj_ratio), adj_ratio_low = mean(adj_ratio_low, adj_ratio_hi = mean(adj_ratio_hi))) %>%
    mutate(cnv_id2 = paste0(sample_id2,'_', cnv_id ))
  
  #count distribuation of CNV pseudotime(adj.ratio)
  ggplot(mf_merged_tbl, aes(x=adj_ratio, fill=deadbody))+
    geom_histogram()+
    scale_fill_manual(values = db_pal)
  
  #age-converted CNVtime per AMP
  db_order <- meta_dt %>% select(deadbody, age) %>% unique() %>% arrange(age) %>% pull(deadbody)
  mf_merged_tbl$deadbod <- factor(mf_merged_tbl$deadbody, levels = db_order)
  cnv_id2_order <- mf_merged_tbl %>% arrange(match(deadbody, db_order), CNV_time)%>% pull(cnv_id2)
  tmp_tbl <- mf_merged_tbl %>% group_by(deadbody,age) %>% count() %>% arrange(match(deadbody, db_order))
  tmp_tbl$cum_sum <- cumsum(tmp_tbl$n)
  tmp_tbl <- tmp_tbl %>% mutate(start_pos = cum_sum -n +1, end_pos=cum_sum)
  nrow(mf_merged_tbl)
  ggplot(mf_merged_tbl, aes(x=cnv_id2, y=CNV_time, color=deadbody))+
    geom_point()+
    geom_errorbar(aes(ymin = CNV_time_low, ymax = CNV_time_hi))+
    geom_segment(data=tmp_tbl, aes(x=start_pos, xend = end_pos, y=age, yend=age), linetype="dashed")+
    scale_x_discrete(limits= cnv_id2_order)+
    scale_color_manual(values = db_pal)+
    coord_flip()
  
  mf_merged_tbl %>% mutate(dist = age - CNV_time) %>% filter(dist >= 10)
  mf_merged_tbl %>% filter(CNV_time < 30)
}


#draw VAF CN plot examples
if(F){
  #early example, 5_LT-4_DS-1_001B10, chr7, CN 2-1-1 dist=159.1Mbps,
  seg_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_seg/5_LT-4_DS-1_001B10.seqz_segments_clean.txt', col_types = cols(`#CHROM`='c'))
  seg_dt <- seg_dt %>% filter(`#CHROM` %in% 6:8)
  seg_dt <- seg_dt %>% mutate(adj.majCN = ifelse(majCN == minCN, majCN+0.1, majCN), adj.minCN = ifelse(majCN == minCN, minCN-0.1, minCN))
  pt_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/5_LT-4_DS-1_001B10.edited.vcf', comment= '##', col_types = cols(`#CHROM`='c'))
  pt_dt <- pt_dt %>% rename('RD_AD_AF' = `5_LT-4_DS-1_001B10`) %>% separate(RD_AD_AF, c('RD','AD','VAF'), sep=':', convert = T) %>% mutate(VAF = VAF*0.01) %>%
    filter(`#CHROM` %in% 6:8)
  pt_dt <- pt_dt %>% mutate(h_line1 = ifelse(`#CHROM` %in% c('6','8'), 0.5,0.33), h_line2 = ifelse(`#CHROM` %in% c('6','8'), NA, 0.66))
  g1 <- ggplot(seg_dt)+
    geom_segment(aes(x=start_pos, xend=end_pos,  y=adj.majCN, yend=adj.majCN), color = "red", size=3)+
    geom_segment(aes(x=start_pos, xend=end_pos, y=adj.minCN, yend=adj.minCN), color = "blue", size=3)+
    facet_grid(cols = vars(`#CHROM`), scales = 'free_x', space = 'free_x')+
    ylab("Copy number")+
    scale_y_continuous(limits = c(0,4))+
    theme_syp+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  g2 <- ggplot(pt_dt, aes(x=POS, y=VAF))+
    geom_point(size=2, alpha=0.5)+
    geom_hline(aes(yintercept = h_line1), linetype ='dashed')+
    geom_hline(aes(yintercept = h_line2), linetype = 'dashed')+
    facet_grid(cols = vars(`#CHROM`), scales = 'free_x', space = 'free_x')+
    scale_y_continuous(limit=c(0,1))+
    theme_syp+theme(axis.text.x = element_text(angle=45, hjust=1))
  plot_grid(g1, g2, ncol = 1, align="v", rel_heights = c(1,1.5))
  
  #late example, 8_ARLL2-3_001E8, chrX, CN2-1-1, dist 155.2Mbps
  seg_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_seg/8_ARLL2-3_001E8.seqz_segments_clean.txt', col_types = cols(`#CHROM`='c'))
  seg_dt <- seg_dt %>% filter(`#CHROM` %in% c(19:22,'X'))
  seg_dt <- seg_dt %>% mutate(adj.majCN = ifelse(majCN == minCN, majCN+0.1, majCN), adj.minCN = ifelse(majCN == minCN, minCN-0.1, minCN))
  pt_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/8_ARLL2-3_001E8.edited.vcf', comment= '##', col_types = cols(`#CHROM`='c'))
  pt_dt <- pt_dt %>% rename('RD_AD_AF' = `8_ARLL2-3_001E8`) %>% separate(RD_AD_AF, c('RD','AD','VAF'), sep=':', convert = T) %>% mutate(VAF = VAF*0.01) %>%
    filter(`#CHROM` %in% c(19:22,'X'))
  pt_dt <- pt_dt %>% mutate(h_line1 = ifelse(`#CHROM` %in% 19:22, 0.5,0.33), h_line2 = ifelse(`#CHROM` %in% 19:22, NA, 0.66))
  g1 <- ggplot(seg_dt)+
    geom_segment(aes(x=start_pos, xend=end_pos,  y=adj.majCN, yend=adj.majCN), color = "red", size=3)+
    geom_segment(aes(x=start_pos, xend=end_pos, y=adj.minCN, yend=adj.minCN), color = "blue", size=3)+
    facet_grid(cols = vars(`#CHROM`), scales = 'free_x', space = 'free_x')+
    ylab("Copy number")+
    scale_y_continuous(limits = c(0,4))+
    theme_syp+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  g2 <- ggplot(pt_dt, aes(x=POS, y=VAF))+
    geom_point(size=2, alpha=0.5)+
    geom_hline(aes(yintercept = h_line1), linetype ='dashed')+
    geom_hline(aes(yintercept = h_line2), linetype = 'dashed')+
    facet_grid(cols = vars(`#CHROM`), scales = 'free_x', space = 'free_x')+
    scale_y_continuous(limit=c(0,1))+
    theme_syp+theme(axis.text.x = element_text(angle=45, hjust=1))
  plot_grid(g1, g2, ncol = 1, align="v", rel_heights = c(1,1.5))
}

#draw mut CN plot examples
if(F){
  #early example, 5_LT-4_DS-1_001B10, chr7, CN 2-1-1 dist=159.1Mbps,
  seg_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_seg/5_LT-4_DS-1_001B10.seqz_segments_clean.txt', col_types = cols(`#CHROM`='c'))
  seg_dt <- seg_dt %>% filter(`#CHROM` %in% 6:8)
  seg_dt <- seg_dt %>% mutate(adj.majCN = ifelse(majCN == minCN, majCN+0.1, majCN), adj.minCN = ifelse(majCN == minCN, minCN-0.1, minCN))
  TCF = 1
  pt_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/5_LT-4_DS-1_001B10.edited.vcf.anv.seqzcn', comment= '##', col_types = cols(`#CHROM`='c'))
  pt_dt <- pt_dt %>% rename('RD_AD_AF' = `5_LT-4_DS-1_001B10`) %>% separate(RD_AD_AF, c('RD','AD','VAF'), sep=':', convert = T) %>% mutate(VAF = VAF*0.01) %>%
    mutate(mutCN = VAF*totCN*(totCN*TCF + (1-TCF)*2)/totCN/TCF) %>%
    mutate(adj.majCN = ifelse(majCN == minCN, majCN+0.03, majCN), adj.minCN = ifelse(majCN == minCN, minCN-0.03, minCN)) %>%
    filter(`#CHROM` %in% 6:8)
  pt_dt <- pt_dt %>% mutate(h_line1 = ifelse(`#CHROM` %in% c('6','8'), 0.5,0.33), h_line2 = ifelse(`#CHROM` %in% c('6','8'), NA, 0.66))
  ggplot(pt_dt, aes(x=POS, y=mutCN))+
    geom_hline(aes(yintercept = adj.majCN), color = 'red',size=3)+
    geom_hline(aes(yintercept = adj.minCN), color = 'blue',size=3)+
    geom_point(size=3, alpha=0.7)+
    scale_y_continuous(limits = c(0,4))+
    facet_grid(cols = vars(`#CHROM`), scales = 'free_x', space = 'free_x')+
    theme_syp+theme(axis.text.x = element_text(angle=45, hjust=1))
  
  
  #late example, 8_ARLL2-3_001E8, chrX, CN2-1-1, dist 155.2Mbps
  seg_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_seg/8_ARLL2-3_001E8.seqz_segments_clean.txt', col_types = cols(`#CHROM`='c'))
  seg_dt <- seg_dt %>% filter(`#CHROM` %in% c(19:22,'X'))
  seg_dt <- seg_dt %>% mutate(adj.majCN = ifelse(majCN == minCN, majCN+0.1, majCN), adj.minCN = ifelse(majCN == minCN, minCN-0.1, minCN))
  TCF=1
  pt_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/8_ARLL2-3_001E8.edited.vcf.anv.seqzcn', comment= '##', col_types = cols(`#CHROM`='c'))
  pt_dt <- pt_dt %>% rename('RD_AD_AF' = `8_ARLL2-3_001E8`) %>% separate(RD_AD_AF, c('RD','AD','VAF'), sep=':', convert = T) %>% mutate(VAF = VAF*0.01) %>%
    mutate(mutCN = VAF*totCN*(totCN*TCF + (1-TCF)*2)/totCN/TCF) %>%
    mutate(adj.majCN = ifelse(majCN == minCN, majCN+0.03, majCN), adj.minCN = ifelse(majCN == minCN, minCN-0.03, minCN)) %>%
    filter(`#CHROM` %in% c(19:22,'X'))
  pt_dt <- pt_dt %>% mutate(h_line1 = ifelse(`#CHROM` %in% 19:22, 0.5,0.33), h_line2 = ifelse(`#CHROM` %in% 19:22, NA, 0.66))
  ggplot(pt_dt, aes(x=POS, y=mutCN))+
    geom_hline(aes(yintercept = adj.majCN), color = 'red',size=3)+
    geom_hline(aes(yintercept = adj.minCN), color = 'blue',size=3)+
    geom_point(size=3, alpha=0.7)+
    scale_y_continuous(limits = c(0,4))+
    facet_grid(cols = vars(`#CHROM`), scales = 'free_x', space = 'free_x')+
    theme_syp+theme(axis.text.x = element_text(angle=45, hjust=1))
}


#Example of Kataegis (EDFig 5a)


source('/home/users/sypark/01_Python_files/circos_R/CIRCOS_R_SYPark.R')
t_id = '2_RT3Fb3_001G8'
print(t_id)
t_gender_ori = meta_dt$gender[meta_dt$sample_id == t_id]
t_gender = gsub('F','XX',gsub('M','XY',t_gender_ori))
t_gender2 = gsub('M','male',gsub('F','female',t_gender_ori))

SNV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/',t_id,'.edited.vcf'), comment = '##', col_types = cols(`#CHROM` = 'c')) %>%
  rename(CHR = `#CHROM`)
SNV_dt <- SNV_dt %>% separate(t_id, c('RD','AD','vaf'), sep=':', convert = T) %>% mutate(vaf = vaf*0.01)
#serial SNV -> DNV
m_SNV_dt <- SNV_dt %>% group_by(CHR) %>% mutate(dist = c(NA, diff(POS)))
m_SNV_dt$dist[is.na(m_SNV_dt$dist) == T] <-m_SNV_dt$POS[is.na(m_SNV_dt$dist) == T]
m_SNV_dt <- m_SNV_dt %>% mutate(var_id = paste(CHR, POS, REF, ALT, sep='_'))
DNV_list <- m_SNV_dt %>% filter(dist == 1) %>% pull(var_id)
for (v_id in DNV_list){
  print(v_id)
  i = which(m_SNV_dt$var_id == v_id)
  m_SNV_dt[i-1,"REF"] <- paste0(m_SNV_dt[i-1,"REF"], m_SNV_dt[i,"REF"])
  m_SNV_dt[i-1,"ALT"] <- paste0(m_SNV_dt[i-1,"ALT"], m_SNV_dt[i,"ALT"])
  m_SNV_dt <- m_SNV_dt[-i,]
}

SV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/structural_variations/SV_final_links/',t_id,'.sv_final.txt'), col_types = cols(`#CHR1`='c', POS1 ='n', CHR2='c', POS2='n')) %>%
  rename(chr1 = `#CHR1`, pos1 = POS1, chr2 = CHR2, pos2 = POS2, ori=terinfo) %>% filter(chr2 != '.')

load(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_extract_RData/',t_id,'.',t_gender,'.100_20_sequenza_extract.RData'))
seqz_dt <- get(gsub('RMCS-34','RMCS_34',paste0(t_id,'.',t_gender,'.100_20_sequenza_extract')))
rm(list=ls(pattern=t_id))



SNV=m_SNV_dt # CHR POS REF ALT vaf color(optional) shape(optional(default pch=16)) size(optional(default cex=0.5))
CNV=CNV_dt # chromosome ("chr1"...) start.pos end.pos A (major CN) B(minor CN) (sequenza segements.txt)
StV=SV_dt # chr1 pos1 chr2 pos2 ori("5to5""3to3""3to5""5to3")
seqz=seqz_dt   # sequenza_extract.RData
purity=1.0
ploidy=2.0
spcs = 'hg19' # hg19 or mm10
gender=t_gender2   # female or male
SNV_plottype = 'rainfall' #c('vaf''rainfall''both')
SNV_height=0.3
CNV_height=0.3
dot_transparency = 0.4
line_transparency = 0.4
min_CNV_scale = 4
seqzR_version = 'v3' # 'v3' or 'v2'
  
pdf(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/structural_variations/',t_id,'for_kataegis_CIRCOS.pdf'),width=12, height=8)

if(T){
  line_pal = c("5to5"=add_transparency("#F23A29", line_transparency), #red
               "3to3"=add_transparency("#F2B705", line_transparency),#orange
               "3to5"=add_transparency("#15A0EC", line_transparency),#blue
               "5to3"=add_transparency("#39A649", line_transparency),#green
               "TRA"=add_transparency("#DD4DF0", line_transparency),
               "Unknown_Mate"=add_transparency("gray", line_transparency))   #purple
  snp_pal = c("C>A"=add_transparency("#15A0EC", dot_transparency), #blue
              "C>G"=add_transparency("#0D0C1B", dot_transparency),#black
              "C>T"=add_transparency("#F23A29", dot_transparency),#red
              "T>A"=add_transparency("#A1A1A1", dot_transparency),#grey
              "T>C"=add_transparency("#5AB440", dot_transparency),#green
              "T>G"=add_transparency("#F2BBC5", dot_transparency),#pink
              'Multi-nucleotide variant'=add_transparency("#DD4DF0", dot_transparency),#purple
              'Indel' = add_transparency("yellow", dot_transparency)) #yellow
  
  if (spcs == 'hg19'){
    chrlist=paste0("chr",c(1:22))
  } else if (spcs == 'mm10'){
    chrlist = paste0("chr", c(1:19))
  }
  if (gender == 'male'){
    chrlist=c(chrlist,"chrX","chrY")
  } else{
    chrlist=c(chrlist,"chrX")
  }
  
  #Draw circos
  circos.clear()
  circos.par("start.degree" = 90, # start from 12 oclock direction
             cell.padding   =c(0.00, 1.00, 0.00, 1.00), # exact ylim will be applied
             track.margin = c(0.012, 0.012), # increased gap between track
             gap.degree = c(rep(1,length(chrlist)-1),4),
             canvas.xlim = c(-1.2,1), canvas.ylim=c(-1, 1)) #left margin for legend
  # cell.padding   =c(0.02, 1.00, 0.02, 1.00)) # original margin
  # Ideogram
  # circos.initializeWithIdeogram()
  circos.initializeWithIdeogram(plotType = NULL, chromosome.index = chrlist, species = spcs) 
  #circos.info()
  circos.genomicIdeogram(track.height = 0.05)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    chr = get.cell.meta.data("sector.index") %>% str_replace("chr","")
    xcenter = get.cell.meta.data("xcenter")
    ycenter = get.cell.meta.data("ylim")[2]
    circos.text(xcenter, ycenter, chr,adj = c(0.5,-1.0)) # adj = position of  chromosome label
    # circos.genomicAxis(h = "top") # for ticks and axis (coordinates)
  })
  
  
  # SNV tract (CHR POS REF ALT vaf)
  SNV$chr = gsub("^", "chr",gsub("chr","",SNV$CHR))
  SNV <- SNV %>% mutate(type = ifelse(nchar(REF) > 1 | nchar(ALT) >1, 
                                      ifelse(nchar(REF) == nchar(ALT), 'Multi-nucleotide variant', 'Indel'),
                                      mgsub::mgsub(paste0(substr(REF,1,1), '>',substr(ALT,1,1)),c("G>T","G>C","G>A","A>T","A>G","A>C"),c("C>A","C>G","C>T","T>A","T>C","T>G"))
  ))
  SNV_dist_v=vector()
  SNV <- SNV %>% arrange(CHR, POS)
  for (n in 1:nrow(SNV)){
    if (n == 1 | as.character(SNV[n,"CHR"]) != as.character(SNV[n-1,"CHR"])){
      SNV_dist_v = c(SNV_dist_v, c(as.numeric(SNV[n,"POS"])))
    } else {
      SNV_dist_v = c(SNV_dist_v, as.numeric(SNV[n,"POS"])- as.numeric(SNV[n-1,"POS"]))
    }
  }
  SNV$dist = log10(SNV_dist_v)
  if (('color' %in% colnames(SNV))== F){
    SNV$color = snp_pal[SNV$type]
  }
  if (('shape' %in% colnames(SNV))==F){
    SNV$shape = 16
  }
  if (('size' %in% colnames(SNV)) == F){
    SNV$size = 0.5
  }
  SNV <- SNV %>% filter(chr %in% chrlist)
  
  if (SNV_plottype == 'rainfall'){
    Draw_dist_plot(SNV,SNV_height,chrlist)
  } else if (SNV_plottype == 'vaf'){
    Draw_vaf_plot(SNV,SNV_height,chrlist)
  } else if (SNV_plottype == 'both'){
    Draw_vaf_plot(SNV,SNV_height/2,chrlist)
    Draw_dist_plot(SNV,SNV_height/2,chrlist)
  }
  
  
  cnv_ymax = ceiling(max(CNV$A[CNV$end.pos-CNV$start.pos > 1E6]))+1
  cnv_ymax = max(min_CNV_scale, cnv_ymax)
  
  # SV and CNV tract
  circos.track(factors=chrlist,
               ylim=c(0,cnv_ymax),
               track.height = CNV_height,
               bg.border="gray")
  StV$chr1 = gsub("^", "chr",gsub("chr","",StV$chr1))
  StV$chr2 = gsub("^", "chr",gsub("chr","",StV$chr2))
  StV <- StV %>% filter(chr1 %in% chrlist & chr2 %in% c(chrlist,'chr.'))
  if(nrow(StV) >0){
    for (i in 1:nrow(StV)){
      if(StV$chr1[i] == StV$chr2[i]){
        # h=ifelse(StV$ori[i] %in% c("5'-5'","3'-3'"), height1,height2)
        Col=line_pal[StV$ori[i]]
      } else if (StV$chr2[i] == 'chr.'){
        Col=line_pal['Unknown_Mate']
      }else {
        # h=NULL
        Col=line_pal["TRA"]
      }
      circos.segments(sector.index = StV$chr1[i],
                      x0 = StV$pos1[i], y0 = 0,
                      x1 = StV$pos1[i], y1 = cnv_ymax,
                      col = Col,
                      lwd=0.7)
      circos.segments(sector.index = StV$chr2[i],
                      x0 = StV$pos2[i], y0 = 0,
                      x1 = StV$pos2[i], y1 = cnv_ymax,
                      col = Col,
                      lwd=0.7)
    }
  }
  
  if (seqzR_version == 'v3'){
    cnt <- function(dr, pu, pl, adr, cnn) {
      (dr/adr*(pu*pl/2+1-pu) - (1-pu))/pu*cnn
    }
  }else if (seqzR_version == 'v2'){
    cnt <- function(dr, pu,pl, adr){
      (dr/adr - (1-pu))*pl/pu
    }
  }
  
  
  if (seqzR_version == 'v3'){
    for(CHR in chrlist) {
      {
        for(lvl in 0:cnv_ymax){ 
          circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="grey", sector.index=CHR)
        }
        seqz.window <- seqz$ratio[[str_replace(CHR,"chr","")]]
        seqz.window <- seqz.window[seqz.window$N >= 1, ]
        seqz.window$mean <- cnt(seqz.window$mean, purity, ploidy, seqz$avg.depth.ratio,2)
        seqz.window$q0 <- cnt(seqz.window$q0, purity, ploidy, seqz$avg.depth.ratio,2)
        seqz.window$q1 <- cnt(seqz.window$q1, purity, ploidy, seqz$avg.depth.ratio,2)
      }
      seqz.window$q1 <- pmin(seqz.window$q1,cnv_ymax)
      seqz.window$mean <- pmin(seqz.window$mean,cnv_ymax)
      seqz.window$q0 <- pmax(seqz.window$q0,0)
      circos.segments(y0 = seqz.window$mean,y1 = seqz.window$mean, x0 = seqz.window$start,
                      x1 = seqz.window$end, lty = 1, lwd = 1, col = "black",sector.index=CHR)
    }
  } else if (seqzR_version == 'v2'){
    for(CHR in chrlist) {
      {
        for(lvl in 0:cnv_ymax){ 
          circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="grey", sector.index=CHR)
        }
        seqz.window <- seqz$ratio[[str_replace(CHR,"chr","")]]
        seqz.window <- seqz.window[seqz.window$N >= 1, ]
        seqz.window$mean <- cnt(seqz.window$mean, purity, ploidy, mean(seqz$gc$adj[,2]))
        seqz.window$q0 <- cnt(seqz.window$q0, purity, ploidy, mean(seqz$gc$adj[,2]))
        seqz.window$q1 <- cnt(seqz.window$q1, purity, ploidy, mean(seqz$gc$adj[,2]))
      }
      seqz.window$q1 <- pmin(seqz.window$q1,cnv_ymax)
      seqz.window$mean <- pmin(seqz.window$mean,cnv_ymax)
      seqz.window$q0 <- pmax(seqz.window$q0,0)
      circos.segments(y0 = seqz.window$mean,y1 = seqz.window$mean, x0 = seqz.window$start,
                      x1 = seqz.window$end, lty = 1, lwd = 1, col = "black",sector.index=CHR)
    }
  }
  
  
  circos.yaxis(side="left",at=(0:cnv_ymax),sector.index="chr1",labels.cex=0.5,tick.length = 0.1)
  
  baseline1=get.cell.meta.data("cell.bottom.radius")
  baseline2=get.cell.meta.data("cell.bottom.radius")
  height1=0.1
  height2=0.1
  if(nrow(StV)>0){
    for (i in 1:nrow(StV)){
      if(StV$chr1[i] == StV$chr2[i]){
        h=ifelse(StV$ori[i] %in% c("5to5","3to3"), height1,height2)
        Col=line_pal[StV$ori[i]]
      } else {
        h=NULL
        Col=line_pal["TRA"]
      }
      circos.link(sector.index1 = StV$chr1[i], StV$pos1[i], 
                  sector.index2 = StV$chr2[i], StV$pos2[i],
                  col = Col,
                  h=h,
                  lwd=0.7,
                  rou=ifelse(StV$ori[i] %in% c("5to5","3to3"), baseline1,baseline2))
    }
  }
  #legend
  # point mutation
  lgd_points = Legend(at = names(snp_pal), type = "points", 
                      legend_gp = gpar(col = snp_pal), title_position = "topleft", 
                      title = paste0("Point mutations\n(n = ", nrow(SNV),")"))

  # SV link
  lgd_sv = Legend(at = c("Head-to-head inversion-type","Tail-to-tail inversion-type","Deletion-type","Duplication-type"), type = "lines", 
                  legend_gp = gpar(col = line_pal[1:4], lwd = 2), title_position = "topleft", 
                  title = paste0("Structural variations\n(n = ",nrow(StV),")"))
  
  lgd_list_vertical = packLegend(lgd_points, lgd_sv)
  draw(lgd_list_vertical, x = unit(0.05, "npc"), y = unit(0, "npc"), just = c("left", "bottom"))
}

dev.off()

