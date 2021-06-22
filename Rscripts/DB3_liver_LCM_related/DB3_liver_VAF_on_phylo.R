# liver LCM data on phylo tree
library(tidyverse)


#load dataset
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210227.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210227.txt')
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


#define functions
Color_VAF_on_phylo_3 <- function(input_tbl, nwk_path, title,legend){  #required columns in in input_tbl: lineage_id, rank, n_vars_same_rank, vaf, n_pointmt
  #color setting
  total_break=100
  start_value=0
  end_value=50
  change_value1=5
  change_value2=20
  col_breaks = seq(start_value, end_value, length.out=total_break)
  a <- length(col_breaks[col_breaks < change_value1] )  #change color 
  b <- length(col_breaks[col_breaks < change_value2])-a  #change color
  rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
  rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
  rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
  my_palette<- c(rampCol1, rampCol2, rampCol3)
  #print legend
  if(legend =='on'){
    if(T){
      plot.new()
      plot(c(0,2),c(start_value,end_value),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, start_value, 0.3,change_value1)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, change_value1, 0.3,change_value2)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, change_value2, 0.3,end_value)
      text(x=0.5, y = c(start_value, change_value1, change_value2, end_value), labels =c(start_value, change_value1, change_value2,  end_value), cex=2)
    }
  }
  #make vaf tbl
  input_tbl$vaf <- round(as.numeric(input_tbl$vaf),2)
  m_input_tbl <- input_tbl %>% mutate(var_g_id = paste0(lineage_id,'_',rank)) %>% filter(is.na(rank)==F) %>% group_by(lineage_id, rank, n_vars_same_rank) %>%
    summarise(mean_vaf = mean(vaf))
  vaf_split_tbl <- tibble(lineage_id = as.character(), vaf_rank = as.character(), mean_vaf=as.numeric())
  for (i in 1:nrow(m_input_tbl)){
    t_id <- m_input_tbl$lineage_id[i]
    t_rank <- m_input_tbl$rank[i]
    t_same <- m_input_tbl$n_vars_same_rank[i]
    t_vaf <- m_input_tbl$mean_vaf[i]
    vaf_ids <- paste0('VAF',seq(t_rank, t_rank + t_same -1))
    t_tbl <- tibble(lineage_id = t_id, vaf_rank = vaf_ids, mean_vaf = t_vaf)
    t_tbl$mean_vaf <- as.numeric(t_tbl$mean_vaf)
    vaf_split_tbl <- bind_rows(vaf_split_tbl, t_tbl)
  }
  annot_tbl <- vaf_split_tbl %>% spread(key=vaf_rank, value=mean_vaf)
  current_cols <- annot_tbl %>% select(starts_with("VAF")) %>% colnames()
  current_cols
  additional_cols <- setdiff(paste0("VAF", 1:20), current_cols)
  if(length(additional_cols) >0){
    annot_tbl <- cbind(annot_tbl, tibble(additional_cols) %>% mutate(value=NA) %>% spread(key=additional_cols, value=value)) %>% as.tibble()
  }
  annot_tbl <- annot_tbl %>% select(lineage_id, paste0('VAF', 1:20))
  annot_tbl <- left_join(annot_tbl, input_tbl %>% select(lineage_id, n_pointmt) %>% unique())
  annot_tbl <- annot_tbl %>% mutate(lid =lineage_id)
  #draw colored phylogenetic tree
  library(ggtree)
  tree <- read.tree(nwk_path)
  g <- ggtree(tree) %<+% annot_tbl + theme_tree2()+
    coord_cartesian(xlim = c(-1,30))+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    geom_text2(aes(x=x-n_pointmt, y = y+0.5, label = lid), hjust=0,color="red",size=2.5)+
    scale_color_gradientn(limits=c(start_value,end_value),colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))+
    ggtitle(title)
  g2 <- ggtree(tree) %<+% annot_tbl+ theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=2)+
    geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    geom_text2(aes(x=x-n_pointmt, y = y+0.5, label = lid), hjust=0,color="red",size=2.5)+
    scale_color_gradientn(limits = c(start_value, end_value), colors = my_palette, breaks = col_breaks)+
    scale_x_continuous(breaks = seq(0,30,1))+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))+
    ggtitle(title)
  print(g2)
}

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


#load data!!!!!!
if(T){
  f_dt <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/db3_liver_WGS/08_DB3_previous_call_readc/DB3_late_liverLCM_f_dt.rds')
  #R10 merge
  #merge R10 data
  f_dt$DB3_liver_R10_dp <- f_dt %>% select(starts_with('DB3_liver_R10')) %>% select(ends_with('dp')) %>% rowSums()
  f_dt$DB3_liver_R10_var <- f_dt %>% select(starts_with('DB3_liver_R10')) %>% select(ends_with('var')) %>% rowSums()
  f_dt$DB3_liver_R10_ref <- f_dt %>% select(starts_with('DB3_liver_R10')) %>% select(ends_with('ref')) %>% rowSums()
  f_dt$DB3_liver_R10_vafpct <- f_dt$DB3_liver_R10_var*100/f_dt$DB3_liver_R10_dp
  f_dt <- f_dt %>% select(-contains('lane'))
  
  #load early data
  input_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3_liver_WGS/08_DB3_previous_call_readc/DB3_early_pointmts.txt.25sCall')
  nrow(input_tbl)
  input_tbl <- input_tbl %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) 
  input_tbl <- input_tbl %>% mutate(lineage_id = gsub('DB3_','',lineage_id2))
  #merge R10 data
  input_tbl$DB3_liver_R10_dp <- input_tbl %>% select(starts_with('DB3_liver_R10')) %>% select(ends_with('dp')) %>% rowSums()
  input_tbl$DB3_liver_R10_var <- input_tbl %>% select(starts_with('DB3_liver_R10')) %>% select(ends_with('var')) %>% rowSums()
  input_tbl$DB3_liver_R10_ref <- input_tbl %>% select(starts_with('DB3_liver_R10')) %>% select(ends_with('ref')) %>% rowSums()
  input_tbl$DB3_liver_R10_vafpct <- input_tbl$DB3_liver_R10_var*100/input_tbl$DB3_liver_R10_dp
  input_tbl <- input_tbl %>% select(-contains('lane'))
}



#Find major lineage in early
t_id = 'C4'
tmp<- input_tbl %>% select(var_id, contains(paste0(t_id,'_')), lineage_id2, -contains('ref'), -contains('_var'), -contains('_dp')) 
colnames(tmp) <- c('var_id', 'vafpct','lineage_id2')
tmp %>% arrange(desc(vafpct)) %>% View()
tmp <- tmp %>% group_by(lineage_id2) %>% summarise(mean_vafpct = mean(vafpct), max_vafpct = max(vafpct)) 
tmp <- left_join(tmp, lng_dt %>% mutate(expected_VAF = early_desc_n*50/trunk_n) %>% select(lineage_id2, expected_VAF))
tmp <- tmp %>% mutate(mean_vaf_ratio = mean_vafpct/expected_VAF, max_vaf_ratio = max_vafpct/expected_VAF)
tmp %>% arrange(desc(max_vaf_ratio))
#tmp %>% arrange(desc(mean_vaf_ratio))

#Find major lineage in late
t_id = 'C3'
tmp<- f_dt %>% select(var_id, contains(paste0(t_id,'_')), lineage_id2, -contains('ref'), -contains('_dp')) 
tmp
colnames(tmp) <- c('var_id', 'var','vafpct','lineage_id2')
tmp %>% arrange(desc(var), desc(vafpct))
tmp <- tmp %>% group_by(lineage_id2) %>% summarise(mean_vafpct = mean(vafpct), max_vafpct = max(vafpct)) 
tmp <- left_join(tmp, lng_dt %>% mutate(expected_VAF = early_desc_n*50/trunk_n) %>% select(lineage_id2, expected_VAF))
tmp <- tmp %>% mutate(mean_vaf_ratio = mean_vafpct/expected_VAF, max_vaf_ratio = max_vafpct/expected_VAF)
tmp %>% arrange(desc(max_vaf_ratio))
#tmp %>% arrange(desc(mean_vaf_ratio))


#plot each sample on the tree
id_list <- gsub('_dp','',input_tbl %>% select(ends_with('dp')) %>% colnames())
id_list <- gsub('DB3_liver_','',id_list)
#pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/DB3_liver_LCM_related/DB3_liver_phylo_total.pdf')
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/DB3_liver_LCM_related/DB3_liver_phylo_var2m.pdf')
n=0
for (t_id in id_list){
  n=n+1
  print(t_id)
  
  ff_dt <- f_dt %>% filter(get(paste0('DB3_liver_',t_id,'_var')) >=2) #filter n_var <2 in late data
  #ff_dt <- f_dt #total
  tmp_dt <- ff_dt %>% select(var_id, lineage_id, contains(paste0(t_id,'_')), -ends_with('ref'), -ends_with('var'),-ends_with('dp'))
  colnames(tmp_dt) <- c('var_id', 'lineage_id','vaf')
  tmp_dt <- tmp_dt %>% group_by(lineage_id) %>% mutate(rank = rank(-vaf, ties.method='first'))
  tmp_dt$n_vars_same_rank =1
  tmp_dt <- tmp_dt %>% filter(rank <20)
  
  s_input_tbl <- input_tbl %>% select(var_id, lineage_id, contains(paste0(t_id,'_')), -ends_with('ref'), -ends_with('var'), -ends_with('dp'))
  colnames(s_input_tbl) <- c('var_id','lineage_id','vaf')
  s_input_tbl <- s_input_tbl %>% group_by(lineage_id) %>% mutate(rank = rank(-vaf, ties.method='first'))
  s_input_tbl <- s_input_tbl %>% mutate(n_vars_same_rank = 1)
  
  s_input_tbl <- bind_rows(tmp_dt,s_input_tbl)
  s_input_tbl <- left_join(s_input_tbl, lng_dt %>% filter(deadbody == 'DB3') %>% select(lineage_id, n_pointmt)) #add n_pointmt info from lng_dt
  nwk_path <- nwk_tbl$nwk_path[nwk_tbl$deadbody=='DB3']
  if (n==1){
    Color_VAF_on_phylo_3(s_input_tbl, nwk_path,paste0('DB3_liver_',t_id), legend='on')
  } else{
    Color_VAF_on_phylo_3(s_input_tbl, nwk_path,paste0('DB3_liver_',t_id), legend='off')
  }
}
dev.off()


#plotting average VAF on phylogenetic tree (main figure 3)
tmp1 <- f_dt %>% select(var_id, lineage_id, ends_with('vafpct'), -starts_with('blood'))
tmp2 <- input_tbl %>% select(var_id, lineage_id, ends_with('vafpct'), -starts_with('blood'))
vaf_dt <- bind_rows(tmp1,tmp2)
tmp3<- vaf_dt %>% select(ends_with('vafpct')) %>% mutate(mean_vafpct = rowMeans(.))
tmp4 <- vaf_dt %>% select(ends_with('vafpct')) %>% mutate(max_vafpct = do.call(pmax, (.))) 
vaf_dt <- do.call(bind_cols, c(vaf_dt %>% select(var_id, lineage_id), tmp3 %>% select(mean_vafpct), tmp4 %>% select(max_vafpct)))
m_vaf_dt <- vaf_dt %>% group_by(lineage_id) %>% mutate(rank = rank(-mean_vafpct, ties.method='first'))
m_vaf_dt$n_vars_same_rank = 1
m_vaf_dt <- left_join(m_vaf_dt, lng_dt %>% filter(deadbody == 'DB3') %>% select(lineage_id, n_pointmt, time_of_lineage))
m_vaf_dt <- m_vaf_dt %>% dplyr::rename(vaf=mean_vafpct)
ggplot(m_vaf_dt)+
  geom_histogram(aes(x=vaf))+
  scale_x_continuous(limits=c(0,5))
fm_vaf_dt <- m_vaf_dt %>% filter(time_of_lineage == 'early' | vaf >=3)
nwk_path <- nwk_tbl$nwk_path[nwk_tbl$deadbody=='DB3']
#pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/DB3_liver_LCM_related/DB3_liver_22s_meanVAF.pdf')
Color_VAF_on_phylo_3(fm_vaf_dt, nwk_path,title='DB3 liver 22sample average', legend='on')
#dev.off()

mm_dt <- vaf_dt %>% group_by(lineage_id) %>% summarise(mean_vafpct = mean(mean_vafpct, na.rm=T), max_vafpct= max(max_vafpct, na.rm=T)) 
mm_dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/DB3_liver_LCM_related/Mean_Max_VAF_perLineage.tsv')

#Draw legend
#color setting
total_break=100
start_value=0
end_value=50
change_value1=5
change_value2=20
col_breaks = seq(start_value, end_value, length.out=total_break)
a <- length(col_breaks[col_breaks < change_value1] )  #change color 
b <- length(col_breaks[col_breaks < change_value2])-a  #change color
rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
my_palette<- c(rampCol1, rampCol2, rampCol3)
#print legend
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/DB3_liver_LCM_related/meanVAF_legend.pdf')
if(T){
    plot.new()
    plot(c(0,2),c(start_value,end_value),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
    colfunc1 <- rev(rampCol1)
    legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
    rasterImage(legend_image1, 0, start_value, 0.3,change_value1)
    colfunc2 <- rev(rampCol2)
    legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
    rasterImage(legend_image2, 0, change_value1, 0.3,change_value2)
    colfunc3 <- rev(rampCol3)
    legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
    rasterImage(legend_image3, 0, change_value2, 0.3,end_value)
    text(x=0.5, y = c(start_value, change_value1, change_value2, end_value), labels =c(start_value, change_value1, change_value2,  end_value), cex=2)
  }
dev.off()


