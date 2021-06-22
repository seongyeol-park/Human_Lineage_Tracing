# DB2 and DB5 HConly filter
library(ComplexHeatmap)
library(circlize)
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db2/04_haplotypecaller/DB2_28s_merged_subs.HC.txt.sampleN.loh_id.g_fi.sam5.txt.28sCall')
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db5/05_haplotypecaller/DB5_28s_merged_subs.HC.txt.sampleN.loh_id.g_fi.sam5.txt.28sCall')
colnames(dt)
dt <- dt %>% mutate(snv_id = paste0(`#CHROM`,':',POS,REF,'>',ALT))
dt <- dt %>% mutate(Nr1=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,2]), Nr3=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,4]))
mVAF <- dt %>% select(ends_with('vafpct')) %>% apply(., 1, function(x) mean(x[x>0 & is.na(x)==F]))
dt$meanVAF <- mVAF
nNA <- dt %>% select(ends_with('vafpct')) %>% apply(.,1,function(x) sum(is.na(x)))
dt$nNA <- nNA

f_dt <- dt %>% filter(Nr1-Nr3 < 10 & Nr3 -casNo < 10 & nNA < 10 & Nr3 >=5 & Nr1 < 28 & mVAF >=15)

hm_mx <- f_dt %>% select(snv_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
colnames(hm_mx) <- gsub('_vafpct','',colnames(hm_mx))
dim(hm_mx)
Heatmap(hm_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")))


