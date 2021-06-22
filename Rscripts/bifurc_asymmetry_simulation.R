library(tidyverse)

lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200722.txt')

bifurc_tbl <- lng_dt %>% filter(branch_n == 2 & early_desc_n > 4)
bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-1')][1]),
                                    sub_desc_n2 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-2')][1]))
bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L1'][1]), sub_desc_n1),
                                    sub_desc_n2 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L2'][1]), sub_desc_n2)) 
bifurc_tbl <- bifurc_tbl %>% mutate(bifurc_ratio = apply(bifurc_tbl[c("sub_desc_n1","sub_desc_n2")], 1, function(x) max(x[1], x[2])/min(x[1], x[2])))
bifurc_tbl <- bifurc_tbl %>% mutate(ratio_fill = ifelse(bifurc_ratio > 7, 7, bifurc_ratio))


sim_res <- list()
sim_res
for (n in 1:100000){
  if(n %% 10000 == 0){
    print(n)
  }
  sim_res[[n]] <- vector()
  for (i in 1:nrow(bifurc_tbl)){
    tot_n <- as.numeric(bifurc_tbl[i,"early_desc_n"])
    n1 <- rbinom(1,tot_n,0.5)
    this_r <- max(tot_n - n1, n1)/min(tot_n - n1, n1)
    sim_res[[n]] <- c(sim_res[[n]], this_r )
  }
}
sim_merged_v <- Reduce(c,sim_res)
sim_median_v <- Reduce(c,lapply(sim_res, median)) 
saveRDS(sim_res, file='/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Asymmetry_simulation/sim_res.1e5.rds')
saveRDS(sim_median_v, file='/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Asymmetry_simulation/sim_median_v.1e5.rds')
