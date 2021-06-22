##Draw_tsne_per_deadbody

this_db = 'DB10'
fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP > 100 & tissue_class != 'cancer_tissue' & time_of_lineage=='early' & deadbody == this_db)
n_total <- fm_vaf_tbl$sample_id %>% unique() %>% length()
available_ids <- fm_vaf_tbl %>% group_by(var_id) %>% dplyr::count() %>% arrange(n) %>% filter(n > n_total*0.5) %>% pull(var_id)
fm_vaf_tbl <- fm_vaf_tbl %>% filter(var_id %in% available_ids)

cutoff1=5;cutoff2=0.5

#3group
diff_tbl <- fm_vaf_tbl  %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early',
                                                                                                                             ifelse(medVAF < cutoff2*0.01, 'late', 'mid')))

#2group
diff_tbl <- fm_vaf_tbl %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early','late'))

Draw_phylo_GroupByVAF(fm_vaf_tbl,
                      lng_dt,
                      this_db,
                      5,
                      0.5)

Draw_Heatmap_Diff(fm_vaf_tbl,
                  diff_tbl)

Draw_tSNE_Diff(fm_vaf_tbl,
               diff_tbl, 
               color_type="AC")

Draw_tSNE_total(fm_vaf_tbl,
                color_type='AC')




if(F){ #3group comparison
  Draw_phylo_timing_TotalAlign_3group(fm_vaf_tbl,lng_dt, this_db, cutoff1, cutoff2)
  diff_tbl_3group <- fm_vaf_tbl  %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early',
                                                                                                                                      ifelse(medVAF < cutoff2*0.01, 'late', 'mid')))
  Draw_Heatmap_Diff_3group(fm_vaf_tbl, diff_tbl_3group)
  Draw_tSNE_Diff_3group(fm_vaf_tbl, diff_tbl_3group)
  
}
if(F){ #2group comparison
  Draw_phylo_timing_TotalAlign_2group(fm_vaf_tbl, lng_dt, this_db, cutoff1)
  diff_tbl_2group <- fm_vaf_tbl %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early','late'))
  Draw_Heatmap_Diff_2group(fm_vaf_tbl, diff_tbl_2group)
  Draw_tSNE_Diff_2group(fm_vaf_tbl, diff_tbl_2group)
  
}