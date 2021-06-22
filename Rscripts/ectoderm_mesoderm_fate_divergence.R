stop()

library(tidyverse)
library(ggsci)
library(scales)
library(ggtree)
library(ggpubr)
library(cowplot)

#load data
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/merged_m_vaf_tbl.rds')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt') %>% filter(current_final_for_lineage == 'Y')

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


#color setting
pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
show_col(pal_combi[c(20,22,23,24)])
group_pal <- pal_combi[c(20,22,12, 23,24)]
names(group_pal) <- c('mixed','fresh_frozen_epidermis','formalin_fixed_epidermis','endoderm','mesoderm')
dl2_pal <- pal_combi[c(8, 12, 4, 3, 13)]
show_col(dl2_pal)
names(dl2_pal) <- c('ectoderm','ecto_mesoderm','mesoderm','meso_endoderm','endoderm')
db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)
lr_pal = c(pal_combi[8], pal_combi[4])
names(lr_pal) = c('rt','lt')
cc_pal = c(pal_combi[24], pal_combi[23])
names(cc_pal) = c('cranio','caudal')
#theme setting
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), strip.text = element_text(size=15), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))

#################################################germ layer###################################################

#2layer plot
if(F){
  this_db = 'DB10'
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group = ifelse(dominant_layer != 'ectoderm', dominant_layer, paste0(tissue_type, '_epidermis'))) 
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(group %in% c('mesoderm','fresh_frozen_epidermis'))
  med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
  med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
  med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
  var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
  #log scale
  ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.001), color = group))+
    geom_point(aes(shape = first_group), size=2, alpha=0.5)+
    geom_errorbar(aes(ymin=log10(q0+0.001), ymax=log10(q1+0.001)), alpha=0.5)+
    geom_line(aes(group = group), alpha=0.5)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    theme_syp+theme(axis.text.x = element_blank())
  #log scale L1 and L2 respectively
  ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.001), color = group))+
    geom_point(aes(shape = first_group), size=2, alpha=0.5)+
    geom_errorbar(aes(ymin=log10(q0+0.001), ymax=log10(q1+0.001)), alpha=0.5)+
    geom_line(aes(group = group), alpha=0.5)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    theme_syp+theme(axis.text.x = element_blank())+
    facet_wrap(~ first_group)
  #normal scale
  ggplot(med_tbl, aes(x=var_id, y=medVAF, color = group))+
    geom_point(aes(shape = first_group), size=2, alpha=0.5)+
    geom_errorbar(aes(ymin=q0, ymax=q1), alpha=0.5)+
    geom_line(aes(group = group), alpha=0.5)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    theme_syp+theme(axis.text.x = element_blank())
  #mesoderm - fresh_frozen_epidermi VAF distance
  m_med_tbl <- med_tbl %>% select(var_id, step_n,group, medVAF) %>% spread(key=group, value=medVAF) %>% mutate(dist = mesoderm - fresh_frozen_epidermis)
  ggplot(m_med_tbl, aes(x=var_id, y=dist))+
    geom_point( size=2, alpha=0.5)+
    scale_x_discrete(limits = var_order)+
    #scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    theme(axis.text.x = element_text(angle=45, hjust=1))
}

# 3group plot
if(F){
  this_db = 'DB9'
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = ifelse(tissue_class2 %in% c('lung', 'liver'), 'endoderm', 
                                                      ifelse(dominant_layer != 'ectoderm', dominant_layer, paste0(tissue_type, '_epidermis'))))
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm', 'fresh_frozen_epidermis'))
  med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
  med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
  med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
  var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
  ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.001), color = group2))+
    geom_point(size=2, alpha=0.5)+
    geom_errorbar(aes(ymin=log10(q0+0.001), ymax=log10(q1+0.001)), alpha=0.5)+
    geom_line(aes(group = group2), alpha=0.5)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
}

# dominant layer2
if(T){
  this_db = 'DB10'
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer2)
  #fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','mixed','fresh_frozen_epidermis'))
  ct_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% count() %>% group_by(group2) %>% summarise(mean_n = mean(n))
  print(ct_tbl)
  incl_ids <- ct_tbl %>% filter(mean_n >=2) %>% pull(group2)
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% incl_ids)
  med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
  med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
  med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
  var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
  ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.001), color = group2))+
    geom_point(size=3, alpha=0.5)+
    stat_smooth(aes(group = group2))+
    coord_cartesian(ylim=log10(c(0,0.5)+0.001))+
    scale_color_manual(values=dl2_pal)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(expand = c(0,0.05), labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
  ggplot(med_tbl, aes(x=var_id, y=medVAF, color = group2))+
    geom_point(size=3, alpha=0.5)+
    stat_smooth(aes(group = group2))+
    scale_color_manual(values=dl2_pal)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(limits=c(0,0.05))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
}

#Proportion of Ectoderm and MesoEndoderm
this_db = 'DB9'
fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer2)
ct_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% count() %>% group_by(group2) %>% summarise(mean_n = mean(n))
print(ct_tbl)
incl_ids <- ct_tbl %>% filter(mean_n >=2) %>% pull(group2)
fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('ectoderm','mesoderm','meso_endoderm','endoderm'))
med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
tmp <- med_tbl %>% select(var_id, group2, medVAF) %>% spread(key=group2, value=medVAF) %>% mutate(MesoEctoratio = ectoderm/mesoderm)
f_tmp <- tmp %>% filter(ectoderm >= 0.003 & mesoderm >= 0.003)
ggplot(f_tmp, aes(MesoEctoratio))+
  geom_histogram(binwidth=0.02)+
  geom_density()+
  ggtitle(this_db)+
  xlab("Mesoderm - Ectoderm VAF ratio")+
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,2))


# x axis with mutation number
if(F){
  this_db = 'DB10'
  f_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP >= 100) 
  f_tbl <- f_tbl %>% mutate(dominant_layer2 = ifelse(tissue_class2 %in% c("lung","liver"), 'endoderm',dominant_layer))%>% filter(dominant_layer != 'ectoderm' | tissue_type == 'fresh_frozen') %>%
    filter(dominant_layer2 != 'mixed')
  med_tbl <- f_tbl  %>% group_by(lineage_id2, start_mtn, rank, n_vars_same_rank, dominant_layer2) %>% summarise(medVAF = median(VAF)) %>%
    filter(dominant_layer2 != 'mixed') 
  ggplot(med_tbl)+
    geom_segment(aes(x=start_mtn+rank-1, xend=start_mtn+rank-1+n_vars_same_rank, y=log10(medVAF+0.001), yend=log10(medVAF+0.001), color = dominant_layer2), size=2, alpha=0.7)+
    scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    ggtitle(t_id)
  m_med_tbl <- med_tbl %>% ungroup() %>% spread(key=dominant_layer2, value=medVAF) %>% mutate(ecto_meso_ratio = ectoderm/mesoderm)
  ggplot(m_med_tbl, aes(x=start_mtn +rank -1, y=ecto_meso_ratio))+
    geom_point()+
    coord_cartesian(ylim=c(0,1))
  ggplot(f_tbl, aes(x=assigned_mean_mtn, y=log10(VAF+0.001), color=dominant_layer2))+
    geom_point(alpha=0.1)+
    geom_smooth()+
    scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))
  el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
  pdf(paste0('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/',this_db,'_VAF_perGroup_perLineage.pdf'))
  for(i in 1:length(el_ids)){
    t_id = el_ids[i]
    print(t_id)
    tmp_v <- unlist(strsplit(t_id, '-'))
    t_id_list <- c()
    for (i in 1:length(tmp_v)){
      t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
    }
    med_tbl <- f_tbl %>% filter(lineage_id2 %in% t_id_list) %>% group_by(lineage_id2, start_mtn, rank, n_vars_same_rank, dominant_layer2) %>% 
      summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4]) %>%
      filter(dominant_layer2 != 'mixed')
    
    g <- ggplot(med_tbl, aes(color=dominant_layer2, fill = dominant_layer2))+
      geom_segment(aes(x=start_mtn+rank-1, xend=start_mtn+rank-1+n_vars_same_rank, y=log10(medVAF+0.001), yend=log10(medVAF+0.001)), size=2, alpha=0.7)+
      geom_rect(aes(xmin=start_mtn+rank-1, xmax=start_mtn+rank-1+n_vars_same_rank, ymin=log10(q0+0.001), ymax=log10(q1+0.001)), alpha=0.3, color=NA)+
      scale_y_continuous(limits = log10(c(0,0.5)+0.001),labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
      scale_x_continuous(limits=c(0,25))+
      ggtitle(t_id)
    print(g)
  }
  dev.off()
}


#wilcox test per individual variant
if(F){
  this_db = 'DB6'
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group = ifelse(dominant_layer != 'ectoderm', dominant_layer, paste0(tissue_type, '_epidermis'))) 
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = ifelse(tissue_class2 %in% c('lung', 'liver'), 'endoderm', group))
  f_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','mixed','fresh_frozen_epidermis'))
  res_tbl = tibble(var_id = as.character(), wilcox_p =as.numeric())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_id)){
    n=n+1
    a <- f_vaf_tbl %>% filter(var_id == t_varid & group == 'fresh_frozen_epidermis') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_id == t_varid & group != 'fresh_frozen_epidermis') %>% pull(VAF)
    if (length(a) > 0){
      res <- wilcox.test(a,b, alternative = 'less')
      pval <- res$p.value
    } else {
      pval <- NA
    }
    res_tbl[n,] <- list(t_varid, pval)
  }
  #distribution of wilcox pvalue
  ggplot(res_tbl)+
    geom_histogram(aes(x=wilcox_p), binwidth=0.01)+
    geom_vline(xintercept = 0.05, linetype = 'dashed')
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% select(lineage_id2, var_id, rank_ctie_first) %>% unique())
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -rank_ctie_first) %>% spread(key=rank2, value=wilcox_p)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  #order in median_vaf
  if(F){
    m_res_tbl <- left_join(fm_vaf_tbl %>% filter(dominant_layer == 'mesoderm') %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(rank_VAF = rank(medVAF, ties.method ='first')), res_tbl)
    m_res_tbl <- m_res_tbl %>% ungroup() %>% select(-var_id, -medVAF) %>% mutate(rank_VAF = paste0('VAF', rank_VAF)) %>% spread(key=rank_VAF, value = wilcox_p)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c("red", "orange"))(n=a)
    rampCol2 <- colorRampPalette(c("orange","yellow"))(n=b)
    rampCol3 <- colorRampPalette(c("yellow", "gray"))(n=total_break-(a+b))
    my_palette<- c(rampCol1, rampCol2, rampCol3)
    
    #legend..................................
    plot.new()
    plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
    colfunc1 <- rampCol1
    legend_image1 <- as.raster(matrix(colfunc1(5), ncol=1))
    rasterImage(legend_image1, 0, 0, 0.3,5)
    colfunc2 <- rampCol2
    legend_image2 <- as.raster(matrix(colfunc2(5), ncol=1))
    rasterImage(legend_image2, 0, 5, 0.3,10)
    colfunc3 <- rampCol3
    legend_image3 <- as.raster(matrix(colfunc3(40), ncol=1))
    rasterImage(legend_image3, 0, 10, 0.3,50)
    text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
    
    
  }
  if(F){ #discrete variable
    cutoff_v = c(0.05)
    change_to_group <- function(numb,cutoff_v){
      if(is.na(numb)){
        return(NA)
      }else{
        n_cutoff_v <- sort(c(numb,cutoff_v), decreasing = T)
        return(paste0('G',min(which(n_cutoff_v == numb))))
      }
    }
    m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), list(~ map_chr(., ~ change_to_group(.x, cutoff_v)))) 
  }
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    #geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    #geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    #geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    #geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    #geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    #geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    #geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    #geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    #geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    #geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    #geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    #geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    #geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    #geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    #geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    #geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    #geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    #ggtitle(paste0(this_db, ', median VAF > ',cutoff,'%, both derm-epi, only derm, none'))+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  }



#wilcox test per variant group
pval_co=0.05
this_db='DB6'
Wilcox_Test_on_Phylo_Ecto <- function(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co){
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group = ifelse(dominant_layer != 'ectoderm', dominant_layer, paste0(tissue_type, '_epidermis'))) 
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = ifelse(tissue_class2 %in% c('lung', 'liver'), 'endoderm', group))
  f_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','fresh_frozen_epidermis') & is.na(rank)== F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_Eclow =as.numeric(), wilcox_p_Echi = as.numeric(), wilcox_p_both = as.numeric(),
                   ecto_mean = as.numeric(), mesoendo_mean = as.numeric(), ecto_n = as.integer(), mesoendo_n = as.integer())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & group == 'fresh_frozen_epidermis') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & group != 'fresh_frozen_epidermis') %>% pull(VAF)
    if (length(a) > 0 & length(b) > 0 ){
      res1 <- wilcox.test(a,b, alternative = 'less')
      pval1 <- res1$p.value
      res2 <- wilcox.test(a,b, alternative = 'greater')
      pval2 <- res2$p.value
      res3 <- wilcox.test(a,b, alternative = 'two.sided')
      pval3 <- res3$p.value
      
    } else {
      pval1 <- NA
      pval2 <- NA
      pval3 <- NA
    }
    res_tbl[n,] <- list(t_varid, pval1, pval2, pval3, mean(a), mean(b), length(a), length(b))
  }
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_Eclow < wilcox_p_Echi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  
  #distribution of wilcox pvalue
  ggplot(res_tbl2)+
    geom_histogram(aes(x=wilcox_p), binwidth=0.01)+
    geom_vline(xintercept = 0.05, linetype = 'dashed')
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                             mutate(var_g_id = paste0(lineage_id2, '_', rank))
    )
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  
  
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c("red", "orange"))(n=a)
    rampCol2 <- colorRampPalette(c("orange","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    #legend..................................
    if(F){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm > Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,50)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c(pal_combi[4], pal_combi[2]))(n=a)
    rampCol2 <- colorRampPalette(c(pal_combi[2],"gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(F){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm < Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,50)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  }
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                              mutate(var_g_id = paste0(lineage_id2, '_', rank)))
    
    el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
    changing_lines = tibble()
    for(t_id in el_ids){
      find_v = 'N'
      tmp_v <- unlist(strsplit(t_id, '-'))
      t_id_list <- c()
      for (i in 1:length(tmp_v)){
        t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
      }
      for(n in 1:(length(t_id_list)-1)){
        rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
        if(length(rank_list) == 0){
          next
        }
        for (m in 1:length(rank_list)){
          t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
          if(is.na(t_pv)){
            next
          }
          if (t_pv < pval_co){
            changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
            find_v = 'Y'
            break
          }
        }
        if (find_v == 'Y'){
          break
        }
      }
    }
    changing_lines <- changing_lines %>% unique()
    changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
  }
  my_func <- function(num){
    if (is.na(num)){
      return (NA)
    } else if(num < -49){
      return (-49)
    } else if ( num > 48){
      return (48) 
    } else {
      return (num)
    }
  }
  m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
    rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_text2(aes(subset = is.na(Source3) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = round(VAF1,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = round(VAF2,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = round(VAF3,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = round(VAF4,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    #geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    #geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    #geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    #geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    #geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    #geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    #geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    #geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    #geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    #geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    #geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    #geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    #geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    #geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    #geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    #geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    ggtitle(paste0(this_db, ', Ectoderm-MesoEndoderm divergence'))+
    geom_point(aes(x=point0.05), shape=1)+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  if(F){
    m_res_tbl2 <- left_join(m_res_tbl2, merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F &dominant_layer == 'mesoderm') %>% mutate(var_g_id = paste0(lineage_id2, '_', rank)) %>%
                              group_by(var_g_id) %>% summarise(medVAF = median(VAF)))
    m_res_tbl2 <- m_res_tbl2 %>% mutate(diff_group = ifelse(wilcox_p < 0.05, 'mid', ifelse(medVAF < 0.005, 'late','early'))) %>% filter(is.na(diff_group) == F)
    m_res_tbl2 %>% group_by(diff_group) %>% count()
    Draw_tSNE_Diff_3group_edit(f_vaf_tbl, m_res_tbl2)
  }
  return(changing_lines)
}

#Ectoderm-mesoendoderm bias timing curve!!!
m_ch_dt <- tibble()
for (this_db in c('DB6','DB8','DB9','DB10')){
  changing_lines <- Wilcox_Test_on_Phylo_Ecto(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co)
  m_ch_dt <- bind_rows(m_ch_dt, changing_lines)
}

m_ch_dt
tmp_tbl <- left_join(
  m_ch_dt,
  lng_dt %>% select(deadbody, lineage_id2, early_desc_n)
)
tmp_tbl <- tmp_tbl %>% mutate(new_mtn = start_mtn + rank) %>% group_by(deadbody, new_mtn) %>% summarise(added_lineage_n = sum(early_desc_n))
tmp_tbl <- tmp_tbl %>% mutate(cum_sum = cumsum(added_lineage_n))
tmp_tbl <- bind_rows(tmp_tbl, tibble(deadbody = tmp_tbl$deadbody %>% unique(), new_mtn = 0,added_lineage_n = 0, cum_sum=0))
tmp_tbl <- left_join(tmp_tbl, lng_dt %>% filter(time_of_lineage== 'early_to_late') %>% group_by(deadbody) %>% summarise(total_lineage_n=n()))
ggplot(tmp_tbl, aes(x=new_mtn, y=cum_sum/total_lineage_n, color=deadbody))+
  geom_point(size=3) +
  geom_line(aes(group = deadbody), linetype='dashed')+
  xlab("No. of mutations")+ylab("Proportion of biased lineages")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,35))+
  scale_color_manual(values=db_pal)+
  ggtitle('Ectoderm - Mesoendoderm bias timing')+
  theme_syp

  
#plot distribution of pval 0.05 
m_ch_dt <- m_ch_dt %>% mutate(deadbody = apply(.['lineage_id2'], 1, function(x) unlist(strsplit(x,'_'))[1]))
ggplot(m_ch_dt, aes(x=point0.05, fill=deadbody))+
  geom_histogram()+
  #geom_density(alpha=0.3)+
  
  scale_x_continuous(limits=c(0,35))+
  scale_fill_manual(values = db_pal)+
  theme_syp+
  facet_wrap(~deadbody)

ggplot(m_ch_dt, aes(x=point0.05, fill=deadbody))+
  geom_density(alpha=0.3)+
  scale_x_continuous(limits=c(0,35))+
  scale_fill_manual(values = db_pal)+
  theme_syp

m_ch_dt
m_ch_dt %>% group_by(deadbody) %>% count()
m_ch_dt %>% group_by(deadbody) %>% summarise(med_point = median(point0.05))
rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
rate_dt <- rate_dt %>% rename(deadbody = ID) %>% select(deadbody, type, mean) %>% spread(key=type, value=mean)
tmp_dt <- left_join(m_ch_dt %>% group_by(deadbody) %>% summarise(med_point = median(point0.05)), rate_dt)
tmp_dt %>% mutate(stage = ifelse(med_point < 2*M0, med_point/M0, (med_point - 2*M0)/M1))



#plotting VAF with mutation number
this_db = 'DB6'
fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group = ifelse(dominant_layer != 'ectoderm', dominant_layer, paste0(tissue_type, '_epidermis'))) 
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = ifelse(tissue_class2 %in% c('lung', 'liver'), 'endoderm', group))
f_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','fresh_frozen_epidermis') & is.na(rank)== F)
f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
alpha = 0.0001
f_vaf_tbl %>% as.data.frame() %>% .[1,]
tmp_tbl <- f_vaf_tbl %>% select(sample_id, var_id, DP, var, VAF, logitVAF, lineage_id2, group2, start_mtn, rank, n_vars_same_rank) 
tmp_tbl <- tmp_tbl %>% group_by(sample_id, lineage_id2, rank, start_mtn, n_vars_same_rank, group2) %>% 
  summarise(meanDP = mean(DP), meanVAR = mean(var)) %>% mutate(meanVAF = round(meanVAR/meanDP,4)) %>% ungroup()
tmp_tbl$logitVAF <- log((tmp_tbl$meanVAF+alpha)/(1-tmp_tbl$meanVAF-alpha))
tmp_tbl$logVAF <- log(tmp_tbl$meanVAF + alpha)
m_tbl <- tibble()
for(i in 1:nrow(tmp_tbl)){
  t_rank = tmp_tbl$rank[i]
  t_samen = tmp_tbl$n_vars_same_rank[i]
  for(m in 1:t_samen){
    m_tbl <- bind_rows(m_tbl, tmp_tbl[i,] %>% mutate(new_rank = rank + m -1))
   }
 }
m_tbl <- m_tbl %>% mutate(new_mtn = start_mtn + new_rank, new_group = ifelse(group2 == 'fresh_frozen_epidermis','ectoderm',group2),
                          new_id = paste0(new_group, '_', new_mtn))
x_order <- m_tbl %>% arrange(new_mtn) %>% pull(new_id) %>% unique()
m_tbl[1,] %>% as.data.frame()
ggplot(m_tbl, aes(x=new_id, y=logitVAF, fill = new_group))+
  geom_boxplot()+
  scale_x_discrete(limits = x_order)+
  theme(axis.text.x = element_text(angle=45, hjust=1))

##########################################left-right fate divergence###########################


# group plot
if(F){
  this_db = 'DB9'
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('LE','UE','trunk'), Source_side, 'unknown'))
  
  print(fm_vaf_tbl %>% group_by(Source_side2, var_id) %>% count() %>% group_by(Source_side2) %>% summarise(mean_n= mean(n)))
  med_tbl <- fm_vaf_tbl %>% group_by(deadbody,Source_side2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
  med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
  med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
  var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
  ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.001), color = Source_side2))+
    geom_point(size=2, alpha=0.5)+
    geom_errorbar(aes(ymin=log10(q0+0.001), ymax=log10(q1+0.001)), alpha=0.5)+
    geom_line(aes(group =  Source_side2), alpha=0.5)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
}


pval_co=0.05
s_range='Eonly' #c('EandT','Eonly') #Extremities and trunk or Extremities only
Wilcox_Test_on_Phylo_LR <- function(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range){
  print(this_db)
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  if(s_range == 'EandT'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('LE','UE','trunk'), Source_side, 'unknown'))
  } else if (s_range == 'Eonly'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('LE','UE'), Source_side, 'unknown'))
  } else {
    print('wrong s_range')
    stop()
  }
  print(fm_vaf_tbl %>% group_by(Source_side2, var_id) %>% count() %>% group_by(Source_side2) %>% summarise(mean_n= mean(n)))
  f_vaf_tbl <- fm_vaf_tbl %>% filter(is.na(rank) == F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_rthi =as.numeric(), wilcox_p_lthi = as.numeric(), wilcox_p_both = as.numeric())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'lt') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'rt') %>% pull(VAF)
    if (length(a) > 0 & length(b) > 0 ){
      res1 <- wilcox.test(a,b, alternative = 'less')
      pval1 <- res1$p.value
      res2 <- wilcox.test(a,b, alternative = 'greater')
      pval2 <- res2$p.value
      res3 <- wilcox.test(a,b, alternative = 'two.sided')
      pval3 <- res3$p.value
      
    } else {
      pval1 <- NA
      pval2 <- NA
      pval3 <- NA
    }
    res_tbl[n,] <- list(t_varid, pval1, pval2, pval3)
  }
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_rthi < wilcox_p_lthi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  #res_tbl <- res_tbl %>% mutate(discrete = floor(wilcox_p*100)) 
  #distribution of wilcox pvalue
  ggplot(res_tbl2)+
    geom_histogram(aes(x=wilcox_p), binwidth=0.01)+
    geom_vline(xintercept = 0.05, linetype = 'dashed')
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                             mutate(var_g_id = paste0(lineage_id2, '_', rank))
    )
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  
  
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c("red", "orange"))(n=a)
    rampCol2 <- colorRampPalette(c("orange","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)

    #legend..................................
    if(F){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF > Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,50)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c(pal_combi[4], pal_combi[2]))(n=a)
    rampCol2 <- colorRampPalette(c(pal_combi[2],"gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(F){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF < Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,50)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  }
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                              mutate(var_g_id = paste0(lineage_id2, '_', rank)))
    
    el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
    changing_lines = tibble()
    for(t_id in el_ids){
      find_v = 'N'
      tmp_v <- unlist(strsplit(t_id, '-'))
      t_id_list <- c()
      for (i in 1:length(tmp_v)){
        t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
      }
      for(n in 1:(length(t_id_list)-1)){
        rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
        if(length(rank_list) == 0){
          next
        }
        for (m in 1:length(rank_list)){
          t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
          if(is.na(t_pv)){
            next
          }
          if (t_pv < pval_co){
            changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
            find_v = 'Y'
            break
          }
        }
        if (find_v == 'Y'){
          break
        }
      }
    }
    changing_lines <- changing_lines %>% unique()
    changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
  }
  my_func <- function(num){
    if (is.na(num)){
      return (NA)
    } else if(num < -49){
      return (-49)
    } else if ( num > 48){
      return (48) 
    } else {
      return (num)
    }
  }
  m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
    rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_text2(aes(subset = is.na(Source3) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = round(VAF1,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = round(VAF2,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = round(VAF3,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = round(VAF4,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    #geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    #geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    #geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    #geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    #geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    #geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    #geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    #geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    #geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    #geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    #geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    #geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    #geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    #geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    #geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    #geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    ggtitle(paste0(this_db, ', Left-right axis, ',s_range))+
    geom_point(aes(x=point0.05), shape=1)+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  return(changing_lines)
}

m_ch_dt <- tibble()
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  changing_lines <- Wilcox_Test_on_Phylo_LR(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co,s_range)
  m_ch_dt <- bind_rows(m_ch_dt, changing_lines)
}
m_ch_dt <- m_ch_dt %>% mutate(deadbody = apply(.['lineage_id2'], 1, function(x) unlist(strsplit(x,'_'))[1]))
ggplot(m_ch_dt, aes(x=point0.05, fill=deadbody))+
  geom_density(alpha=0.7)+
  scale_x_continuous(limits=c(0,35))+
  scale_fill_manual(values = db_pal)+
  xlab(paste0('pval < ',pval_co))+
  ggtitle(paste0('Lt-Rt axis discrimination, ', s_range))+
  theme_syp

tmp_tbl <- left_join(
  m_ch_dt,
  lng_dt %>% select(deadbody, lineage_id2, early_desc_n)
)
tmp_tbl <- tmp_tbl %>% mutate(new_mtn = start_mtn + rank) %>% group_by(deadbody, new_mtn) %>% summarise(added_lineage_n = sum(early_desc_n))
tmp_tbl <- tmp_tbl %>% mutate(cum_sum = cumsum(added_lineage_n))
tmp_tbl <- bind_rows(tmp_tbl, tibble(deadbody = tmp_tbl$deadbody %>% unique(), new_mtn = 0,added_lineage_n = 0, cum_sum=0))
tmp_tbl <- left_join(tmp_tbl, lng_dt %>% filter(time_of_lineage== 'early_to_late') %>% group_by(deadbody) %>% summarise(total_lineage_n=n()))
ggplot(tmp_tbl, aes(x=new_mtn, y=cum_sum/total_lineage_n, color=deadbody))+
  geom_jitter(size=3, alpha=0.5, height=0, width=0.5) +
  geom_line(aes(group = deadbody), linetype='dashed', alpha=0.5)+
  xlab("No. of mutations")+ylab("Proportion of biased lineages")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,35))+
  scale_color_manual(values=db_pal)+
  ggtitle('Left-Right bias timing')+
  theme_syp

#Comparing first-occuring mutations (L1_initial vs L2_initial)
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  if(s_range == 'EandT'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('LE','UE','trunk'), Source_side, 'unknown'))
  } else if (s_range == 'Eonly'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('LE','UE'), Source_side, 'unknown'))
  } else {
    print('wrong s_range')
    stop()
  } 
  colnames(fm_vaf_tbl)
  this_lid= 'L1'
  subdt <- fm_vaf_tbl %>% filter(lineage_id2 == paste0(this_db,'_',this_lid) & Source_side2 %in% c('lt','rt') & rank == 1) 
  g1<- ggplot(subdt, aes(x=Source_side2, y=VAF, fill = Source_side2))+
    geom_boxplot()+
    geom_point(size=3, alpha=0.3)+
    ggtitle(paste0(this_db, '_', this_lid))+
    stat_compare_means(comparisons = list(c('lt','rt')), aes(label = ..p.signif..), method="wilcox")+
    scale_fill_manual(values = lr_pal)+
    theme_syp+theme(axis.title.x = element_blank(), legend.position = 'none')
  this_lid= 'L2'
  subdt <- fm_vaf_tbl %>% filter(lineage_id2 == paste0(this_db,'_',this_lid) & Source_side2 %in% c('lt','rt') & rank == 1) 
  g2<- ggplot(subdt, aes(x=Source_side2, y=VAF, fill=Source_side2))+
    geom_boxplot()+
    geom_point(size=3, alpha=0.3)+
    ggtitle(paste0(this_db, '_', this_lid,', ',s_range))+
    stat_compare_means(comparisons = list(c('lt','rt')), aes(label = ..p.signif..),  method="wilcox")+
    scale_fill_manual(values = lr_pal)+
    theme_syp+theme(axis.title.x = element_blank(), legend.position = 'none')
  plot_grid(g1,g2,nrow=1) %>% print()
}
#Number of Lt and Rt clonal samples in initial bifurcation
lng_dt %>% filter(deadbody == 'DB10' & time_of_lineage == 'early_to_late') %>% 
  mutate(Source_side2 = ifelse(Source_class %in% c('LE','UE','trunk') & is.na(Source_side) == F, Source_side, 'unknown'),
         first_group = substr(lineage_id, 1,2)) %>%
  group_by(first_group, Source_side2) %>% count()

##########################################cranio-caudal fate divergence###########################


# group plot
if(F){
  this_db='DB9'
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('UE','HN') | anatomy_class2 == 'chest', 'cranio', 
                                                            ifelse(anatomy_class1 == 'LE' | anatomy_class2 %in% c('abdomen','pubis'), 'caudal','unknown'))) 
  print(fm_vaf_tbl %>% group_by(Source_side2, var_id) %>% count() %>% group_by(Source_side2) %>% summarise(mean_n= mean(n)))
  med_tbl <- fm_vaf_tbl %>% group_by(deadbody,Source_side2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
  med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
  med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
  var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
  ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.001), color = Source_side2))+
    geom_point(size=2, alpha=0.5)+
    geom_errorbar(aes(ymin=log10(q0+0.001), ymax=log10(q1+0.001)), alpha=0.5)+
    geom_line(aes(group =  Source_side2), alpha=0.5)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(labels = c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.001))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
}


pval_co=0.05
s_range = 'HTE' #c('HTE','HandLE') # head trunk extremities, head and LEx
Wilcox_Test_on_Phylo_CC <- function(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range){
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  if(s_range == 'HTE'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('UE','HN') | anatomy_class2 == 'chest', 'cranio',
                                                              ifelse(anatomy_class1 == 'LE' | anatomy_class2 %in% c('abdomen','pubis'), 'caudal','unknown')))
  } else if (s_range == 'HandLE'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('UE','HN') , 'cranio', 
                                                              ifelse(anatomy_class1 == 'LE', 'caudal','unknown')))
  } else {
    print('wrong s_range')
    stop()
  }
 
  print(this_db)
  print(fm_vaf_tbl %>% group_by(Source_side2, var_id) %>% count() %>% group_by(Source_side2) %>% summarise(mean_n= mean(n)))
  f_vaf_tbl <- fm_vaf_tbl %>% filter(is.na(rank) == F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_rthi =as.numeric(), wilcox_p_lthi = as.numeric(), wilcox_p_both = as.numeric())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'caudal') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'cranio') %>% pull(VAF)
    if (length(a) > 0 & length(b) > 0 ){
      res1 <- wilcox.test(a,b, alternative = 'less')
      pval1 <- res1$p.value
      res2 <- wilcox.test(a,b, alternative = 'greater')
      pval2 <- res2$p.value
      res3 <- wilcox.test(a,b, alternative = 'two.sided')
      pval3 <- res3$p.value
      
    } else {
      pval1 <- NA
      pval2 <- NA
      pval3 <- NA
    }
    res_tbl[n,] <- list(t_varid, pval1, pval2, pval3)
  }
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_rthi < wilcox_p_lthi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  #res_tbl <- res_tbl %>% mutate(discrete = floor(wilcox_p*100)) 
  #distribution of wilcox pvalue
  ggplot(res_tbl2)+
    geom_histogram(aes(x=wilcox_p), binwidth=0.01)+
    geom_vline(xintercept = 0.05, linetype = 'dashed')
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                             mutate(var_g_id = paste0(lineage_id2, '_', rank))
    )
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  
  
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c(pal_combi[24], pal_combi[17]))(n=a)
    rampCol2 <- colorRampPalette(c(pal_combi[17],"gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    #legend..................................
    if(F){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Cranio. VAF > Caudal. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,50)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c(pal_combi[23], pal_combi[19]))(n=a)
    rampCol2 <- colorRampPalette(c(pal_combi[19],"gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(F){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF < Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,50)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  }
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                              mutate(var_g_id = paste0(lineage_id2, '_', rank)))
    
    el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
    changing_lines = tibble()
    for(t_id in el_ids){
      find_v = 'N'
      tmp_v <- unlist(strsplit(t_id, '-'))
      t_id_list <- c()
      for (i in 1:length(tmp_v)){
        t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
      }
      for(n in 1:(length(t_id_list)-1)){
        rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
        if(length(rank_list) == 0){
          next
        }
        for (m in 1:length(rank_list)){
          t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
          if(is.na(t_pv)){
            next
          }
          if (t_pv < pval_co){
            changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
            find_v = 'Y'
            break
          }
        }
        if (find_v == 'Y'){
          break
        }
      }
    }
    changing_lines <- changing_lines %>% unique()
    changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
  }
  my_func <- function(num){
    if (is.na(num)){
      return (NA)
    } else if(num < -49){
      return (-49)
    } else if ( num > 48){
      return (48) 
    } else {
      return (num)
    }
  }
  m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
    rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_text2(aes(subset = is.na(Source3) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = round(VAF1,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = round(VAF2,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = round(VAF3,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = round(VAF4,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    #geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    #geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    #geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    #geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    #geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    #geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    #geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    #geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    #geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    #geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    #geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    #geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    #geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    #geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    #geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    #geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    ggtitle(paste0(this_db, ', Cranio-caudal axis, ', s_range))+
    geom_point(aes(x=point0.05), shape=1)+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  return(changing_lines)
}

m_ch_dt <- tibble()
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  changing_lines <- Wilcox_Test_on_Phylo_CC(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co,s_range)
  m_ch_dt <- bind_rows(m_ch_dt, changing_lines)
}
m_ch_dt <- m_ch_dt %>% mutate(deadbody = apply(.['lineage_id2'], 1, function(x) unlist(strsplit(x,'_'))[1]))
ggplot(m_ch_dt, aes(x=point0.05, fill=deadbody))+
  geom_density(alpha=0.7)+
  scale_x_continuous(limits=c(0,35))+
  scale_fill_manual(values = db_pal)+
  xlab(paste0('pval < ',pval_co))+
  theme_syp

tmp_tbl <- left_join(
  m_ch_dt,
  lng_dt %>% select(deadbody, lineage_id2, early_desc_n)
)
tmp_tbl <- tmp_tbl %>% mutate(new_mtn = start_mtn + rank) %>% group_by(deadbody, new_mtn) %>% summarise(added_lineage_n = sum(early_desc_n))
tmp_tbl <- tmp_tbl %>% mutate(cum_sum = cumsum(added_lineage_n))
tmp_tbl <- bind_rows(tmp_tbl, tibble(deadbody = tmp_tbl$deadbody %>% unique(), new_mtn = 0,added_lineage_n = 0, cum_sum=0))
tmp_tbl <- left_join(tmp_tbl, lng_dt %>% filter(time_of_lineage== 'early_to_late') %>% group_by(deadbody) %>% summarise(total_lineage_n=n()))
ggplot(tmp_tbl, aes(x=new_mtn, y=cum_sum/total_lineage_n, color=deadbody))+
  geom_point(size=3) +
  geom_line(aes(group = deadbody), linetype='dashed')+
  xlab("No. of mutations")+ylab("Proportion of biased lineages")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,35))+
  scale_color_manual(values=db_pal)+
  ggtitle('Cranio -Caudal bias timing')+
  theme_syp




#Comparing first-occuring mutations (L1_initial vs L2_initial)
this_db='DB6'
fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('UE','HN') | anatomy_class2 == 'chest', 'cranio', 
                                                          ifelse(anatomy_class1 == 'LE' | anatomy_class2 %in% c('abdomen','pubis'), 'caudal','unknown')))
colnames(fm_vaf_tbl)
this_lid= 'L1'
subdt <- fm_vaf_tbl %>% filter(lineage_id2 == paste0(this_db,'_',this_lid) & Source_side2 %in% c('cranio','caudal') & rank == 1) 
g1<- ggplot(subdt, aes(x=Source_side2, y=VAF, fill = Source_side2))+
  geom_boxplot()+
  geom_point(size=3, alpha=0.3)+
  ggtitle(paste0(this_db, '_', this_lid))+
  stat_compare_means(comparisons = list(c('cranio','caudal')), aes(label = ..p.signif..), method="wilcox")+
  scale_fill_manual(values = cc_pal)+
  theme_syp+theme(axis.title.x = element_blank(), legend.position = 'none')
this_lid= 'L2'
subdt <- fm_vaf_tbl %>% filter(lineage_id2 == paste0(this_db,'_',this_lid) & Source_side2 %in% c('cranio','caudal') & rank == 1) 
g2<- ggplot(subdt, aes(x=Source_side2, y=VAF, fill=Source_side2))+
  geom_boxplot()+
  geom_point(size=3, alpha=0.3)+
  ggtitle(paste0(this_db, '_', this_lid))+
  stat_compare_means(comparisons = list(c('cranio','caudal')), aes(label = ..p.signif..),  method="wilcox")+
  scale_fill_manual(values = cc_pal)+
  theme_syp+theme(axis.title.x = element_blank(), legend.position = 'none')
plot_grid(g1,g2,nrow=1)

#Number of Lt and Rt clonal samples in initial bifurcation
lng_dt %>% filter(deadbody == 'DB10' & time_of_lineage == 'early_to_late') %>% 
  mutate(Source_side2 = ifelse(Source_class %in% c('LE','UE','trunk') & is.na(Source_side) == F, Source_side, 'unknown'),
         first_group = substr(lineage_id, 1,2)) %>%
  group_by(first_group, Source_side2) %>% count()


#########Representative phylogenetic tree of DB6

#1. Ectoderm-MesoEndoderm
this_db = 'DB6'
fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group = ifelse(dominant_layer != 'ectoderm', dominant_layer, paste0(tissue_type, '_epidermis'))) 
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = ifelse(tissue_class2 %in% c('lung', 'liver'), 'endoderm', group))
f_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','fresh_frozen_epidermis') & is.na(rank)== F)
f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
res_tbl = tibble(var_g_id = as.character(), wilcox_p_Eclow =as.numeric(), wilcox_p_Echi = as.numeric(), wilcox_p_both = as.numeric(),
                 ecto_mean = as.numeric(), mesoendo_mean = as.numeric(), ecto_n = as.integer(), mesoendo_n = as.integer())
n=0
for(t_varid in unique(f_vaf_tbl$var_g_id)){
  n=n+1
  a <- f_vaf_tbl %>% filter(var_g_id == t_varid & group == 'fresh_frozen_epidermis') %>% pull(VAF)
  b <- f_vaf_tbl %>% filter(var_g_id == t_varid & group != 'fresh_frozen_epidermis') %>% pull(VAF)
  if (length(a) > 0 & length(b) > 0 ){
    res1 <- wilcox.test(a,b, alternative = 'less')
    pval1 <- res1$p.value
    res2 <- wilcox.test(a,b, alternative = 'greater')
    pval2 <- res2$p.value
    res3 <- wilcox.test(a,b, alternative = 'two.sided')
    pval3 <- res3$p.value
    
  } else {
    pval1 <- NA
    pval2 <- NA
    pval3 <- NA
  }
  res_tbl[n,] <- list(t_varid, pval1, pval2, pval3, mean(a), mean(b), length(a), length(b))
}
res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_Eclow < wilcox_p_Echi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)

#wilcox pvalue on phylogenetic tree
if(T){ #order in rank (random rank in ties)
  m_res_tbl <- left_join(res_tbl,
                         merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                           mutate(var_g_id = paste0(lineage_id2, '_', rank))
  )
  m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
  max_VAF_len =20
  for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
    m_res_tbl <- m_res_tbl %>% mutate(this = NA)
    colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
  }
}

if(T){ #continuous variable
  library(colorRamps)
  total_break=50
  col_breaks = seq(0, 0.5, length.out=total_break)
  a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
  
  b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
  rampCol1 <- colorRampPalette(c("red", "orange"))(n=a)
  rampCol2 <- colorRampPalette(c("orange","gray"))(n=b)
  rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
  my_palette1<- c(rampCol1, rampCol2, rampCol3)
  
  #legend..................................
  if(T){
    plot.new()
    plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm > Ectoderm')
    colfunc1 <- rev(rampCol1)
    legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
    rasterImage(legend_image1, 0, 0, 0.3,5)
    colfunc2 <- rev(rampCol2)
    legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
    rasterImage(legend_image2, 0, 5, 0.3,10)
    colfunc3 <- rev(rampCol3)
    legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
    rasterImage(legend_image3, 0, 10, 0.3,50)
    text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
  }
  
  total_break=50
  col_breaks = seq(0, 0.5, length.out=total_break)
  a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
  
  b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
  rampCol1 <- colorRampPalette(c(pal_combi[4], pal_combi[2]))(n=a)
  rampCol2 <- colorRampPalette(c(pal_combi[2],"gray"))(n=b)
  rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
  my_palette2<- c(rampCol1, rampCol2, rampCol3)
  
  
  #legend..................................
  if(T){
    plot.new()
    plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm < Ectoderm')
    colfunc1 <- rev(rampCol1)
    legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
    rasterImage(legend_image1, 0, 0, 0.3,5)
    colfunc2 <- rev(rampCol2)
    legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
    rasterImage(legend_image2, 0, 5, 0.3,10)
    colfunc3 <- rev(rampCol3)
    legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
    rasterImage(legend_image3, 0, 10, 0.3,50)
    text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.5), cex=2)
  }
  
  my_palette <- c(rev(my_palette2), my_palette1)
  names(my_palette) <- seq(-49,48,1)
}

#assign point changing to <0.05 and >= 0.05
if(T){
  m_res_tbl2 <- left_join(res_tbl2,
                          merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                            mutate(var_g_id = paste0(lineage_id2, '_', rank)))
  
  el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
  changing_lines = tibble()
  for(t_id in el_ids){
    find_v = 'N'
    tmp_v <- unlist(strsplit(t_id, '-'))
    t_id_list <- c()
    for (i in 1:length(tmp_v)){
      t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
    }
    for(n in 1:(length(t_id_list)-1)){
      rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
      if(length(rank_list) == 0){
        next
      }
      for (m in 1:length(rank_list)){
        t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
        if(is.na(t_pv)){
          next
        }
        if (t_pv < pval_co){
          changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
          find_v = 'Y'
          break
        }
      }
      if (find_v == 'Y'){
        break
      }
    }
  }
  changing_lines <- changing_lines %>% unique()
  changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
}
my_func <- function(num){
  if (is.na(num)){
    return (NA)
  } else if(num < -49){
    return (-49)
  } else if ( num > 48){
    return (48) 
  } else {
    return (num)
  }
}
m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
  rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
  coord_cartesian(xlim = c(-1,35))+
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
  ggtitle(paste0(this_db, ', Ectoderm-MesoEndoderm divergence'))+
  scale_color_manual(values = my_palette)+
  theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
print(g)

