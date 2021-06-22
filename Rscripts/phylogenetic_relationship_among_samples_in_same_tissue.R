stop('abc')

# additional analysis for first revision

#libraries

library(tidyverse)
library(ggtree)
library(cowplot)

# load data
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
#DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_200918.txt')
#rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
#pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_201023.txt')



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


plot_list=list()
#lineage distance derived from same tissue (~1cm2)
this_db='DB9'
plot_list[this_db] = list()

dist_dt <- tibble()

for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  print(this_db)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  d = fortify(tree)
  d1 <- d %>% dplyr::rename(x_child = x, y_child=y, node_child=node, node_parent=parent, label_child = label) %>% select(node_parent, node_child, label_child, isTip, x_child, y_child)
  d2 <- d %>% dplyr::rename(node_parent=node, label_parent=label, x_parent=x, y_parent=y) %>% select(node_parent, label_parent, x_parent, y_parent)
  d_merged <- left_join(d1, d2)
  tree_line_size=0.5
  g1 <-ggplot(d_merged)+
    geom_segment(aes(x= x_parent, xend = x_parent, y=y_parent, yend=y_child), size=tree_line_size)+
    geom_segment(aes(x= x_parent, xend = x_child, y=y_child, yend=y_child), size=tree_line_size)+
    coord_cartesian(xlim=c(0,35))+
    theme_bw()+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), title = element_text(size=15, face='bold'))+
    ggtitle(this_db)
  max_v <- max(d_merged$x_child)*1.1
  g2 <- ggplot(d_merged)+
    geom_segment(data = subset(d_merged, isTip == T),aes(x= x_parent, xend = x_parent, y=y_parent, yend=y_child), size=tree_line_size)+
    geom_segment(aes(x= x_parent, xend = x_child, y=y_child, yend=y_child), size=tree_line_size)+
    geom_segment(data = subset(d_merged, isTip == T), aes(x= x_child, xend = max_v, y=y_child, yend=y_child), linetype="dashed", color="gray")+
    coord_cartesian(xlim=c(35,max_v + 100))+
    scale_x_continuous(expand = c(0,0))+
    theme_bw()+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
  
  
  f_meta_dt <- meta_dt %>% filter(deadbody == this_db)
  f_meta_dt$tissue_id %>% unique() %>% length()
  tid_multi <- f_meta_dt %>% group_by(tissue_id) %>% dplyr::count() %>% filter(n>1) %>% pull(tissue_id)
  f_meta_dt <- f_meta_dt %>% filter(tissue_id %in% tid_multi)
  nrow(f_meta_dt)
  f_meta_dt <- left_join(f_meta_dt,d_merged %>% dplyr::rename(lineage_id = label_child) %>% select(lineage_id, y_child))
  f_meta_dt <- left_join(f_meta_dt, lng_dt %>% filter(deadbody == this_db) %>% select(lineage_id, time_of_lineage) %>% unique())
  tmp <- f_meta_dt %>% filter(time_of_lineage == 'late')%>% mutate(mother_lid = lapply(lineage_id, function(x) {
    lid_list = unlist(strsplit(x, '-'))
    paste(lid_list[1:(length(lid_list)-1)], collapse='-')
  }) %>% as.character()) 
  target_mlids <- tmp %>% group_by(tissue_id, mother_lid) %>% count() %>% filter(n>1) %>% pull(mother_lid)
  target_lids <- tmp %>% filter(mother_lid %in% target_mlids) %>% pull(lineage_id)
  
  #tid_both_lineage <- f_meta_dt %>% mutate(lineage_group = substr(lineage_id, 1,2)) %>% group_by(tissue_id, lineage_group) %>% dplyr::count() %>% group_by(tissue_id) %>%
  #  dplyr::count() %>% filter(n>1) %>% pull(tissue_id)
  #f_meta_dt <- f_meta_dt %>% mutate(dist_group = ifelse(tissue_id %in% tid_both_lineage, "both","single"))
  f_meta_dt <- f_meta_dt %>% mutate(dist_group = ifelse(lineage_id %in% target_lids, "late_branched","others"))
  min_max_dt <- f_meta_dt %>% group_by(tissue_id) %>% summarise(min_y=min(y_child), max_y=max(y_child))
  min_max_dt <- left_join(min_max_dt, f_meta_dt %>% select(tissue_id, dist_group))
  
  dist_dt <- bind_rows(dist_dt, f_meta_dt %>% select(sample_id, tissue_id, dist_group))

  
  g3 <- ggplot(f_meta_dt, aes(x=tissue_id, y=y_child, color=dist_group))+
    geom_hline(yintercept = 1:max(d_merged$y_child), linetype='dashed', color='gray')+
    geom_segment(data=min_max_dt, aes(x=tissue_id, xend=tissue_id, y=min_y, yend=max_y), color="black")+
    geom_point(size=3)+
    scale_color_manual(values = c("late_branched"="red", "others"="blue"), labels = c("Late branched", "Othres"))+
    xlab("Small tissues")+
    theme_bw()+
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(),
          panel.border = element_blank(), legend.title = element_blank(), legend.position='none')
  
  plot_list[[this_db]][["tree"]] <- plot_grid(g1, g2, g3, nrow=1, align="h", axis="tb", rel_widths=c(1,1,1.5))
  #plot_list[[this_db]][["tree"]]
}


#count tissues with only late branching clones.
tis_with_late_br <- dist_dt %>% filter(dist_group =="late_branched") %>% pull(tissue_id) 
n_only_late <- dist_dt %>% group_by(tissue_id, dist_group) %>% count() %>% filter(tissue_id %in% tis_with_late_br) %>% group_by(tissue_id) %>% count() %>% filter(n==1) %>% nrow()
print(paste0(n_only_late*100/length(dist_dt$tissue_id %>% unique()),'%'))




#calculate distance score of samples in same tissues
f_meta_dt <- meta_dt %>% filter(deadbody == this_db)
tid_multi <- f_meta_dt %>% group_by(tissue_id) %>% dplyr::count() %>% filter(n>1) %>% pull(tissue_id)
res_tbl <- tibble(tissue_id = as.character(), lineage1 = as.character(), lineage2 = as.character(), shared_mts = as.numeric(), distance_score= as.numeric())
j=0
for (i in 1:length(tid_multi)){
  this_tid <- tid_multi[i]
  lid_list <- f_meta_dt %>% filter(tissue_id == this_tid) %>% pull(lineage_id)
  comb_mx <- combn(lid_list, 2)
  
  for (n in 1:ncol(comb_mx)){
    j = j+1
    print(n)
    str1 <- comb_mx[1,n]
    str1 <- gsub('L','',gsub('-','', str1))
    str2 <- comb_mx[2,n]
    str2 <- gsub('L','',gsub('-','',str2))
    min_length <- min(nchar(str1), nchar(str2))
    cha_list1 <- unlist(strsplit(str1,""))[1:min_length]
    cha_list2 <- unlist(strsplit(str2,""))[1:min_length]
    first_T_index <- min(which(cha_list1 == cha_list2))
    first_F_index <- min(which(cha_list1 != cha_list2))
    
    #caclulate the number of shared mutation
    if(first_T_index > 1){
      shared_mts = 0
    } else{
      shared_id_list <- cha_list1[1:(first_F_index-1)]
      shared_mts=0
      for (m in 1:length(shared_id_list)){
        target_id = paste0("L",paste(shared_id_list[1:m], collapse='-'))
        shared_mts = shared_mts + lng_dt$n_pointmt[lng_dt$lineage_id == target_id & lng_dt$deadbody == this_db]
      } 
    }
    
    distance_score = 1/(shared_mts + 1)
    res_tbl[j,] <- list(this_tid, comb_mx[1,n], comb_mx[2,n], shared_mts, distance_score)
  }
}
total_mean_v <- res_tbl %>% group_by(tissue_id) %>% summarise(mean_score = mean(distance_score)) %>% pull(mean_score) %>% mean()
total_mean_v #0.62

#make random distribution and calculate the distance score
n_random = 1000 #number of performing random selection
r_mean_score_v <- vector()
for (k in 1:n_random){
  if (k %% 100 == 0){
    print(k)
  }
  
  f_meta_dt <- meta_dt %>% filter(deadbody == this_db)
  tid_multi <- f_meta_dt %>% group_by(tissue_id) %>% dplyr::count() %>% filter(n>1) %>% pull(tissue_id)
  ct_tbl <- f_meta_dt %>% filter(tissue_id %in% tid_multi) %>% group_by(tissue_id) %>% dplyr::count()
  total_count <- sum(ct_tbl$n)
  ct_tbl$end_numb <- cumsum(ct_tbl$n)
  ct_tbl <- ct_tbl %>% mutate(start_numb = end_numb - n + 1)
  r_id_list <- sample(f_meta_dt$lineage_id, total_count)  #random selection
  r_tbl <- tibble()
  for (n in 1:nrow(ct_tbl)){
    tissue_id = ct_tbl[n,"tissue_id"] %>% as.character()
    start_numb=ct_tbl[n,"start_numb"] %>% as.numeric()
    end_numb = ct_tbl[n, "end_numb"] %>% as.numeric()
    tmp_tbl <- tibble(lineage_id = r_id_list[start_numb:end_numb], tissue_id = tissue_id)
    r_tbl <- bind_rows(r_tbl, tmp_tbl)
  }

 
  j=0
  res_tbl <- tibble(tissue_id = as.character(), lineage1 = as.character(), lineage2 = as.character(), shared_mts = as.numeric(), distance_score= as.numeric())
  for (i in 1:length(tid_multi)){
    this_tid <- tid_multi[i]
    lid_list <- r_tbl %>% filter(tissue_id == this_tid) %>% pull(lineage_id)
    comb_mx <- combn(lid_list, 2)
    
    for (n in 1:ncol(comb_mx)){
      j = j+1
      str1 <- comb_mx[1,n]
      str1 <- gsub('L','',gsub('-','', str1))
      str2 <- comb_mx[2,n]
      str2 <- gsub('L','',gsub('-','',str2))
      min_length <- min(nchar(str1), nchar(str2))
      cha_list1 <- unlist(strsplit(str1,""))[1:min_length]
      cha_list2 <- unlist(strsplit(str2,""))[1:min_length]
      first_T_index <- min(which(cha_list1 == cha_list2))
      first_F_index <- min(which(cha_list1 != cha_list2))
      
      #caclulate the number of shared mutation
      if(first_T_index > 1){
        shared_mts = 0
      } else{
        shared_id_list <- cha_list1[1:(first_F_index-1)]
        shared_mts=0
        for (m in 1:length(shared_id_list)){
          target_id = paste0("L",paste(shared_id_list[1:m], collapse='-'))
          shared_mts = shared_mts + lng_dt$n_pointmt[lng_dt$lineage_id == target_id & lng_dt$deadbody == this_db]
        } 
      }
      
      distance_score = 1/(shared_mts + 1)
      res_tbl[j,] <- list(this_tid, comb_mx[1,n], comb_mx[2,n], shared_mts, distance_score)
    }
  }
  r_total_mean <- res_tbl %>% group_by(tissue_id) %>% summarise(mean_score = mean(distance_score)) %>% pull(mean_score) %>% mean()
  r_mean_score_v  <- c(r_mean_score_v, r_total_mean)
}
rank_obs <- which(sort(c(r_mean_score_v, total_mean_v)) == total_mean_v)
pval = min(1001-rank_obs, rank_obs)/n_random
score_dt <- tibble(distance_score = r_mean_score_v)
plot_list[[this_db]][["sim"]]  <- ggplot(score_dt, aes(x=distance_score))+
  geom_density(fill="lightgray")+
  scale_x_continuous(limits=c(0,1))+
  geom_vline(xintercept = total_mean_v, color="red")+
  ggtitle(paste0(this_db, ", p-value = ", pval))+
  theme_syp
plot_list[[this_db]][["sim"]]

#saveRDS(plot_list, '/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/Phylogenetic_relationship_same_tissue/plot_tree_simulation.rds')

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig4a.pdf', width=20, height=10)
grid.newpage()
pushViewport(viewport(x=0, y=0.5, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB3']][["tree"]], newpage = F)
popViewport(1)
pushViewport(viewport(x=0.3, y=0.5, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB8']][["tree"]], newpage = F)
popViewport(1)
pushViewport(viewport(x=0.6, y=0, width = 0.4, height = 1, just = c('left', 'bottom')))
print(plot_list[['DB6']][["tree"]], newpage = F)
popViewport(1)
pushViewport(viewport(x=0, y=0, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB9']][["tree"]], newpage = F)
popViewport(1)
pushViewport(viewport(x=0.3, y=0, width = 0.3, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB10']][["tree"]], newpage = F)
popViewport(1)
dev.off()

new_list=list()
n=0
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  n=n+1
  new_list[[n]] <- plot_list[[this_db]][["sim"]]
}
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig4b.pdf', width=20, height=5)
plot_grid(plotlist = new_list, nrow=1, align="h", axis="tb")
dev.off()



##lineage distance derived from same sources
this_db='DB10'
if(T){
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  d = fortify(tree)
  d1 <- d %>% dplyr::rename(x_child = x, y_child=y, node_child=node, node_parent=parent, label_child = label) %>% select(node_parent, node_child, label_child, isTip, x_child, y_child)
  d2 <- d %>% dplyr::rename(node_parent=node, label_parent=label, x_parent=x, y_parent=y) %>% select(node_parent, label_parent, x_parent, y_parent)
  d_merged <- left_join(d1, d2)
  
  g1 <-ggplot(d_merged)+
    geom_segment(aes(x= x_parent, xend = x_parent, y=y_parent, yend=y_child), size=1)+
    geom_segment(aes(x= x_parent, xend = x_child, y=y_child, yend=y_child), size=1)+
    coord_cartesian(xlim=c(0,35))+
    theme_bw()+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), title = element_text(size=15, face='bold'))+
    ggtitle(this_db)
  max_v <- max(d_merged$x_child)*1.1
  g2 <- ggplot(d_merged)+
    geom_segment(aes(x= x_parent, xend = x_parent, y=y_parent, yend=y_child), size=1)+
    geom_segment(aes(x= x_parent, xend = x_child, y=y_child, yend=y_child), size=1)+
    geom_segment(data = subset(d_merged, isTip == T), aes(x= x_child, xend = max_v, y=y_child, yend=y_child), linetype="dashed")+
    coord_cartesian(xlim=c(200,max_v + 100))+
    scale_x_continuous(expand = c(0,0))+
    theme_bw()+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
  
  
  f_meta_dt <- meta_dt %>% filter(deadbody == this_db)
  sid_multi <- f_meta_dt %>% group_by(Source_id) %>% dplyr::count() %>% filter(n>1) %>% pull(Source_id)
  f_meta_dt <- f_meta_dt %>% filter(Source_id %in% sid_multi)
  f_meta_dt <- left_join(f_meta_dt,d_merged %>% dplyr::rename(lineage_id = label_child) %>% select(lineage_id, y_child))
  sid_both_lineage <- f_meta_dt %>% mutate(lineage_group = substr(lineage_id, 1,2)) %>% group_by(Source_id, lineage_group) %>% dplyr::count() %>% 
    group_by(Source_id) %>% dplyr::count() %>% filter(n>1) %>% pull(Source_id)
  f_meta_dt <- f_meta_dt %>% mutate(dist_group = ifelse(Source_id %in% sid_both_lineage, "both","single"))
  min_max_dt <- f_meta_dt %>% group_by(Source_id) %>% summarise(min_y=min(y_child), max_y=max(y_child))
  min_max_dt <- left_join(min_max_dt, f_meta_dt %>% select(Source_id, dist_group))
  
  g3 <- ggplot(f_meta_dt, aes(x=Source_id, y=y_child, color=dist_group))+
    geom_hline(yintercept = 1:max(d_merged$y_child), linetype='dashed')+
    geom_point(size=3)+
    geom_segment(data=min_max_dt, aes(x=Source_id, xend=Source_id, y=min_y, yend=max_y))+
    scale_color_manual(values = c("both"="red", "single"="blue"), labels = c("Both lineage", "Single lineage"))+
    theme_bw()+
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle=45, hjust=1), axis.ticks.y = element_blank(), axis.title = element_blank(),
          panel.border = element_blank(), legend.title = element_blank())
  
  plot_grid(g1, g2, g3, nrow=1, align="h", axis="tb")
}
