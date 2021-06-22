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
  merged_tbl <- bind_rows(merged_tbl, ts_dt %>% filter(!var_id %in% merged_tbl$var_id) %>% mutate(lineage_id2 = paste0(this_db,'_GSNP')))
  return(merged_tbl)
}

Make_VAF_tbl <- function(merged_tbl, this_db){
  vaf_tbl1 <- merged_tbl  #copy table
  #extract dp and var information
  for (c_name in col_list){  
    print(c_name)
    vaf_tbl1[[paste0(c_name, '_dp')]] <- apply(vaf_tbl1[c_name], 1, function(x) as.numeric(gsub('dp=','',unlist(strsplit(x,';'))[1])))
    vaf_tbl1[[paste0(c_name, '_var')]] <- apply(vaf_tbl1[c_name], 1, function(x) as.numeric(gsub('var=','',unlist(strsplit(x,';'))[3])))
    vaf_tbl1 <- vaf_tbl1 %>% select(-c_name)
  }
  
  vaf_tbl1 <- vaf_tbl1 %>% select(var_id, lineage_id2, WGS_bloodVAF, ends_with('_dp'), ends_with('_var'))
  #caclulate binom.confint
  library(binom)
  for (c_name in col_list){ 
    print(c_name)
    vaf_tbl1[[paste0(c_name, '_var_low')]] <- round(vaf_tbl1[[paste0(c_name, '_dp')]]*binom.confint(vaf_tbl1[[paste0(c_name, '_var')]],vaf_tbl1[[paste0(c_name, '_dp')]], method="exact")$lower,4)
    vaf_tbl1[[paste0(c_name, '_var_up')]] <- round(vaf_tbl1[[paste0(c_name, '_dp')]]*binom.confint(vaf_tbl1[[paste0(c_name, '_var')]],vaf_tbl1[[paste0(c_name, '_dp')]], method="exact")$upper,4)
    vaf_tbl1[[paste0(c_name, '_vaf')]] <- round(vaf_tbl1[[paste0(c_name, '_var')]]/vaf_tbl1[[paste0(c_name, '_dp')]],4)
    vaf_tbl1[[paste0(c_name, '_vaf_low')]] <- round(vaf_tbl1[[paste0(c_name, '_var_low')]]/vaf_tbl1[[paste0(c_name, '_dp')]],4)
    vaf_tbl1[[paste0(c_name, '_vaf_up')]] <- round(vaf_tbl1[[paste0(c_name, '_var_up')]]/vaf_tbl1[[paste0(c_name, '_dp')]],4)
  }
  #join with lng_dt
  vaf_tbl1 <- left_join(vaf_tbl1, lng_dt %>% filter(deadbody == this_db) %>% select(lineage_id2, time_of_lineage, early_desc_n,step_n,time_of_lineage, start_mtn, end_mtn))
  tmp_tbl <- vaf_tbl1 %>% select(ends_with('_var'))
  vaf_tbl1$nPositive <- apply(tmp_tbl, 1, function(x) sum(x>0))
  vaf_tbl1$meanVarRC <- round(apply(tmp_tbl, 1, function(x) mean(x[x>0])),0)
  tmp_tbl2 <- vaf_tbl1 %>% select(ends_with('_dp'))
  vaf_tbl1$hidp100 <- apply(tmp_tbl2, 1, function(x) sum(x>=100))
  return(vaf_tbl1)
}

Add_Expected_VAF <- function(vaf_tbl1, this_db){
  #consistency of expected VAF and actual VAF of tissues
  id_list <- vaf_tbl1 %>% filter(is.na(time_of_lineage)==F) %>% .$lineage_id2 %>% unique()
  res_tbl1 = tibble(lineage_id2 = as.character(), expected_VAF = as.numeric())
  n=0
  for (this_lid in id_list){
    n= n+1
    lng_dt %>% filter(time_of_lineage =='early_to_late' & grepl(this_lid, lineage_id2)) %>% nrow() -> a
    lng_dt %>% filter(time_of_lineage =='early_to_late' & deadbody == this_db ) %>% nrow() -> b
    lower_margin = 0.48*binom.confint(a,b, method="exact")$lower
    upper_margin = 0.48*binom.confint(a,b, method="exact")$upper
    expected_value = 0.48*a/b
    res_tbl1[n,] <- list(this_lid, expected_value)
  }
  vaf_tbl1 <- left_join(vaf_tbl1, res_tbl1)
  return(vaf_tbl1)
}
Draw_phylo_GroupByVAF <-function(
  fm_vaf_tbl,
  lng_dt,
  this_db,
  cutoff1=NULL,
  cutoff2=NULL){
  max_VAF_len=20
  early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2 
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(lineage_id2 %in% early_ids)
  rank_tbl <- fm_vaf_tbl %>% group_by(lineage_id2, var_id) %>% summarise(medianVAF = median(VAF)) %>% filter(grepl('germline',lineage_id2)==F) %>% 
    group_by(lineage_id2) %>% mutate(rank = rank(-medianVAF, ties.method = "first")) 
  rank_tbl <- rank_tbl %>% mutate(rank2 = paste0('VAF',rank)) %>% select(-var_id, -rank) %>% spread(key=rank2, value=medianVAF)
  for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(rank_tbl))){
    rank_tbl <- rank_tbl %>% mutate(this = NA)
    colnames(rank_tbl) <- gsub('this',id,colnames(rank_tbl))
  }
  mx <- rank_tbl %>% ungroup() %>% select(starts_with('VAF'))
  mx = mx*100
  if (is.null(cutoff1)){
    mx[mx < cutoff2] <- 0
    mx[mx >= cutoff2] <- 1
  } else{
    mx[mx < cutoff2] <- 0
    mx[mx >= cutoff2 & mx < cutoff1 ] <- 1
    mx[mx >= cutoff1] <- 2
  }
  

  m_rank_tbl <- bind_cols(rank_tbl %>% select(lineage_id2), mx %>% as.tibble())
   m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db),m_rank_tbl)
  m_lng_dt <- m_lng_dt %>% dplyr::rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
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
    ggtitle(paste0(this_db, ', median VAF cutoff1: ',cutoff1,'%, cutoff2: ',cutoff2,'%'))+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}

Draw_phylo_GroupByManual <-function(
  fm_vaf_tbl,
  lng_dt,
  this_db,
  diff_tbl
  ){
  tmp_tbl <- left_join(fm_vaf_tbl, diff_tbl) 
  tmp_tbl <- tmp_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first))
  rank_tbl <- tmp_tbl %>% select(lineage_id2, rank2, diff_group) %>% unique() %>%  spread(key=rank2, value=diff_group)
  max_VAF_len=20
  for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(rank_tbl))){
    rank_tbl <- rank_tbl %>% mutate(this = NA)
    colnames(rank_tbl) <- gsub('this',id,colnames(rank_tbl))
  }
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db),rank_tbl)
  m_lng_dt <- m_lng_dt %>% dplyr::rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
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
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}




Draw_Heatmap_Diff <- function(fm_vaf_tbl,
                              diff_tbl
                              )
  {
  #set colors
  library(ggsci)
  pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
  anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
  names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
  lr_pal = pal_combi[c(8, 11, 10, 18)]
  names(lr_pal) <- c("lt","rt","center","unknown")
  cc_pal = pal_combi[c(24,23,18)]
  names(cc_pal) <- c('cranio','caudal','unknown')

  #make mx
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  e_mx <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  e_mx <- e_mx[rowSums(is.na(e_mx)) < ncol(e_mx)/2,]
  e_mx <- t(scale(t(e_mx)))
  l_mx <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  l_mx <- l_mx[rowSums(is.na(l_mx)) < ncol(l_mx)/2,]
  l_mx <- t(scale(t(l_mx)))
  if("mid" %in% unique(diff_tbl$diff_group)){
    m_mx <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
    m_mx <- m_mx[rowSums(is.na(m_mx)) < ncol(m_mx)/2,]
    m_mx <- t(scale(t(m_mx)))
    mx_list = c('e_mx','m_mx','l_mx')
  }else{
    mx_list = c('e_mx','l_mx')
  }
  name_change = c('early','mid','late')
  names(name_change) <- c('e_mx','m_mx','l_mx')
  #draw heatmap
  annot_dt <- fm_vaf_tbl %>% dplyr::rename(left_right = Source_side, cranio_caudal = Source_side_CC, anatomy = anatomy_class3) %>%
    select(sample_id, left_right, cranio_caudal, anatomy) %>% unique() %>% as.data.frame() %>% column_to_rownames('sample_id')
  anno_colors <- list(left_right = lr_pal, cranio_caudal = cc_pal, anatomy = anatomy_cl3_pal)
  library(pheatmap)
  p_list=list()
  n=1
  for (this_mx in mx_list){
    if(n < length(mx_list)){
      p_list[[n]] <- pheatmap(get(this_mx), annotation_col = annot_dt, show_colnames = F, clustering_method ="ward.D", cluster_rows = T, main=name_change[this_mx], annotation_legend = F, annotation_colors = anno_colors)[[4]]
    } else {
      p_list[[n]] <- pheatmap(get(this_mx), annotation_col = annot_dt, show_colnames = F, clustering_method ="ward.D", cluster_rows = T, main=name_change[this_mx], annotation_legend = T, annotation_colors = anno_colors)[[4]]
    }
    n = n+1
  }
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=1, rel_widths = c(rep(1, length(mx_list)-1),1.25))
}


Draw_PCA_Diff <- function(fm_vaf_tbl, diff_tbl){
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  library(missForest) #NA value imputation using missForest
  Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  Early <- missForest(t(Early))$ximp
  Mid <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
  Mid <- missForest(t(Mid))$ximp
  Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  Late <- missForest(t(Late))$ximp
  p_list=list()
  n=0
  for (timing in c('Early', 'Mid', 'Late')){
    n = n+1
    df_pca <- prcomp(as.data.frame(get(timing)), scale=T, center = T)
    df_out <- as.data.frame(df_pca$x)
    dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), pmeta_dt)
    p_list[[n]] <- ggplot(dt, aes(x=PC1, y=PC2, color=anatomy_class3))+
      geom_point(size=3, alpha=0.7)+
      scale_color_manual(values=anatomy_cl3_pal)+
      theme_syp+
      ggtitle(timing)
  }
  p_list[[1]] <- p_list[[1]]+theme(legend.position = 'none')
  p_list[[2]] <- p_list[[2]]+theme(legend.position = 'none')
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=1, rel_widths = c(1,1,1.25))
  }


Draw_tSNE_Diff <- function(
  fm_vaf_tbl,
  diff_tbl,
  color_type = "LR", #c("LR","CC","AC") LR: left-right, CC: cranio-caudal, AC: anatomy_class
  seed = 1 #the input number of set.seed()
  ){
  set.seed(seed)
  
  #color setting
  library(ggsci)
  pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
  anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
  names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
  lr_pal = pal_combi[c(8, 11, 10, 18)]
  names(lr_pal) <- c("lt","rt","center","unknown")
  cc_pal = pal_combi[c(24,23,18)]
  names(cc_pal) <- c('cranio','caudal','unknown')
  gl_pal = c( "#314da0","#c7ffaf","#4cce75","#409607","#6936af")
  names(gl_pal) <- c("meso_endoderm","mesoderm","ecto_mesoderm","ectoderm","endoderm")
  anatomy_cl1_pal <- c("#c97734","#6e5dd3","#42dd90","#ef7223","#dbb11a")
  names(anatomy_cl1_pal) <- c("internal_organ","HN","trunk","UE","LE")
  anatomy_cl2_pal <- c("#ffdbbc","#86e8d6","#dd58a4","#eeffad","#baf9a4","#e0e241","#3718e2","#dd7573","#21f25c","#b52760","#1845f7","#9693ea","#4e85f4","#f9bdda")
  names(anatomy_cl2_pal) <- c("chest","abdomen","shoulder","lower_arm","leg","pubis","thigh","upper_arm","axilla","elbow","foot","hand","knee","arm")
  #tissue_cl2_pal <- randomColor(21)
  #names(tissue_cl2_pal) <- fm_vaf_tbl$tissue_class2 %>% unique()
  tc3_pal <-c("gray",'#8FDE9E',"#DB9E51","#D9CDE1","#DDE265","#D5938F","#CCD4AD","#B250E1","#83DA94","#83E557","#73B2CF","#E3478B","#81E4D4","#DC89D2")
  names(tc3_pal)<- c("non_internal_organ","liver","colon","large_vessel","heart","stomach","lung","pancreas","intestine","spleen","blood","kidney","ovary","uterus")
  this_pal = c(anatomy_cl3_pal, lr_pal, cc_pal, gl_pal, anatomy_cl1_pal, anatomy_cl2_pal, tc3_pal)

  #make mx in each group
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  mx <- mx[rowSums(mx, na.rm = T)!=0,]
  library(missForest) #NA value imputation using missForest
  Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  Early <- Early[rowSums(is.na(Early)) < ncol(Early)/2,]
  Early <- missForest(t(Early))$ximp
  Early <- Early[rowSums(Early) >0,]
  Early <- t(scale(t(Early)))
  Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  Late <- Late[rowSums(is.na(Late)) < ncol(Late)/2,]
  Late <- missForest(t(Late))$ximp
  Late <- Late[rowSums(Late) >0,]
  Late <- t(scale(t(Late)))
  if("mid" %in% diff_tbl$diff_group %>% unique){
    Mid <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
    Mid <- Mid[rowSums(is.na(Mid)) < ncol(Mid)/2,]
    Mid <- missForest(t(Mid))$ximp
    Mid <- Mid[rowSums(Mid) >0,]
    Mid <- t(scale(t(Mid)))
  }
  
  #run tsne
  library(Rtsne)
  p_list=list()
  n=0
  if("mid" %in% diff_tbl$diff_group %>% unique){
    timing_group = c("Early","Mid","Late")
  } else{
    timing_group = c("Early","Late")
  }
  for (timing in timing_group){
    print(timing)
    n = n+1
    df_tsne <- Rtsne(get(timing), perplexity=10)
    df_out <- as.data.frame(df_tsne$Y)
    rownames(df_out) <- rownames(get(timing))
    colnames(df_out) <- c('tSNE1','tSNE2')
    dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), pmeta_dt)
    if(color_type== "LR"){
      dt <- dt %>% mutate(color_type = Source_side)
    } else if (color_type == 'CC'){
      dt <- dt %>% mutate(color_type = Source_side_CC)
    } else if (color_type == 'AC'){
      dt <- dt %>% mutate(color_type = anatomy_class3)
    } else if (color_type == 'GL'){
      dt <- dt %>% mutate(color_type = dominant_layer2)
    } else if (color_type == 'AC1'){
      dt <- dt %>% mutate(color_type = anatomy_class1)
    } else if (color_type == 'AC2'){
      dt <- dt %>% mutate(color_type = anatomy_class2)
    } else if (color_type == 'TC3'){
      dt <- dt %>% mutate(color_type = tissue_class3)
    }
    if(n < length(timing_group)){
      p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=color_type))+
        geom_point(size=5, alpha=0.7)+
        scale_color_manual(values=this_pal)+
        theme_syp+
        theme(legend.position = 'none')+
        ggtitle(timing)  
    } else{
      p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=color_type))+
        geom_point(size=5, alpha=0.7)+
        scale_color_manual(values=this_pal)+
        theme_syp+
        ggtitle(timing)  
    }
  }
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=1, rel_widths = c(rep(1,length(timing_group)-1),1.25))
}

Draw_tSNE_total <- function(
  fm_vaf_tbl,
  color_type = "LR",  #c("LR","CC","AC") LR: left-right, CC: cranio-caudal, AC: anatomy_class
  seed=1
){
  set.seed(seed)
  
  #color setting
  library(ggsci)
  pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
  anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
  names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
  lr_pal = pal_combi[c(8, 11, 10, 18)]
  names(lr_pal) <- c("lt","rt","center","unknown")
  cc_pal = pal_combi[c(24,23,18)]
  names(cc_pal) <- c('cranio','caudal','unknown')
  this_pal = c(anatomy_cl3_pal, lr_pal, cc_pal)
  #make mx in each group
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  mx <- mx[rowSums(mx, na.rm = T)!=0,]
  mx <- mx[rowSums(is.na(mx)) < ncol(mx)/2,]
  library(missForest) #NA value imputation using missForest
  mx <- missForest(t(mx))$ximp
  mx <- mx[rowSums(mx) >0,]
  mx <- t(scale(t(mx)))
  #run tsne
  library(Rtsne)
  df_tsne <- Rtsne(mx, perplexity=10)
  df_out <- as.data.frame(df_tsne$Y)
  rownames(df_out) <- rownames(mx)
  colnames(df_out) <- c('tSNE1','tSNE2')
  dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), pmeta_dt)
  if(color_type== "LR"){
    dt <- dt %>% mutate(color_type = Source_side)
  } else if (color_type == 'CC'){
    dt <- dt %>% mutate(color_type = Source_side_CC)
  } else if (color_type == 'AC'){
    dt <- dt %>% mutate(color_type = anatomy_class3)
  }
  g<- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=color_type))+
        geom_point(size=5, alpha=0.7)+
        scale_color_manual(values=this_pal)+
        theme_syp
  print(g)
}


Wilcox_Test_on_Phylo_Ecto_MesoEndo <- function(
  f_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, legend='off'
  ){
  print(this_db)
  fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db)
  f_vaf_tbl <- fm_vaf_tbl %>% filter(dominant_layer2 %in% c('mesoderm','endoderm','ectoderm') & is.na(rank)== F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_Eclow =as.numeric(), wilcox_p_Echi = as.numeric(), wilcox_p_both = as.numeric(),
                   ecto_mean = as.numeric(), mesoendo_mean = as.numeric(), ecto_n = as.integer(), mesoendo_n = as.integer())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){  #group id
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 == 'ectoderm') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 != 'ectoderm') %>% pull(VAF)
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
                           f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
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
    rampCol1 <- colorRampPalette(c("#E65996", "#ED8BB6"))(n=a)
    rampCol2 <- colorRampPalette(c("#ED8BB6","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(legend == 'on'){
      plot.new()
      plot(c(0,2),c(0,10),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm > Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      #colfunc3 <- rev(rampCol3)
      #legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      #rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c("#245DB5", "#668ECC"))(n=a)
    rampCol2 <- colorRampPalette(c("#668ECC","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(legend == 'on'){
      plot.new()
      plot(c(0,2),c(0,10),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm < Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      #colfunc3 <- rev(rampCol3)
      #legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      #rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
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
  m_lng_dt <- m_lng_dt %>% dplyr::rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
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
    ggtitle(this_db)+
    #ggtitle(paste0(this_db, ', Ectoderm-MesoEndoderm divergence'))+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15), plot.title = element_text(size=20), plot.margin= unit(c(0.2,1,0.2,0.2),"cm"))
  print(g)
}

Wilcox_Test_on_Phylo_Endo_Meso <- function(f_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co){
  print(this_db)
  fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db)
  f_vaf_tbl <- fm_vaf_tbl %>% filter(dominant_layer2 %in% c('mesoderm','endoderm') & is.na(rank)== F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_Enlow =as.numeric(), wilcox_p_Enhi = as.numeric(), wilcox_p_both = as.numeric(),
                   endo_mean = as.numeric(), meso_mean = as.numeric(), endo_n = as.integer(), meso_n = as.integer())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){  #group id
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 == 'endoderm') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 != 'endoderm') %>% pull(VAF)
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
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_Enlow < wilcox_p_Enhi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
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
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Mesoderm > Ecdoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
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
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Mesoderm < Endoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  }
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
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
  m_lng_dt <- m_lng_dt %>% dplyr::rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
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
    ggtitle(paste0(this_db, ', Endoderm-Mesoderm divergence'))+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  #return(changing_lines)
}

Wilcox_Test_on_Phylo_Ecto_Meso <- function(f_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co){
  print(this_db)
  fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db)
  f_vaf_tbl <- fm_vaf_tbl %>% filter(dominant_layer2 %in% c('mesoderm','ectoderm') & is.na(rank)== F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_Enlow =as.numeric(), wilcox_p_Enhi = as.numeric(), wilcox_p_both = as.numeric(),
                   endo_mean = as.numeric(), meso_mean = as.numeric(), endo_n = as.integer(), meso_n = as.integer())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){  #group id
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 == 'ectoderm') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & dominant_layer2 != 'ectoderm') %>% pull(VAF)
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
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_Enlow < wilcox_p_Enhi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
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
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Mesoderm > Ecdoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
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
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Mesoderm < Endoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
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
  m_lng_dt <- m_lng_dt %>% dplyr::rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
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
    ggtitle(paste0(this_db, ', Ectoderm-Mesoderm divergence'))+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
}


Wilcox_Test_on_Phylo_LR <- function(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range, legend='off'){
  print(this_db)
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db)
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
    rampCol1 <- colorRampPalette(c('#023858', "#3690c0"))(n=a)
    rampCol2 <- colorRampPalette(c("#3690c0","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    #legend..................................
    if(legend == 'on'){
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
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c('#df1919',"#fcd57b"))(n=a)
    rampCol2 <- colorRampPalette(c("#fcd57b","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(legend == 'on'){
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
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
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
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,30))+
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
    ggtitle(this_db)+
    #ggtitle(paste0(this_db, ', Left-right axis, ',s_range))+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
}

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
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_crhi =as.numeric(), wilcox_p_cahi = as.numeric(), wilcox_p_both = as.numeric())
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
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_crhi < wilcox_p_cahi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
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
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF > Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
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
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF < Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
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
    coord_cartesian(xlim = c(-1,30))+
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
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  return(changing_lines)
}
