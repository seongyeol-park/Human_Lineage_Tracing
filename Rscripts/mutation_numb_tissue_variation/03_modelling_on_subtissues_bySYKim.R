##
## Extract only a subset of tissues that have at least two cells
## check whether fitting on these subsets give similar conclusions. 
## 
## Su Yeon
## 7/8/2020 
#######################################################


require(data.table)
require(readxl)
require(RColorBrewer)

library(foreach)
library(doParallel)

require(ggplot2)
require(grid)
require(gridExtra)
require(pheatmap)

require(tidyverse)
require(ggsci)

require(lme4)
require(lmerTest)

require(VCA)


require(randomcoloR)



#########################################################
## Directory setting 
#########################################################


prjName <- "05_lineagetracing"
cat("prjName\t", prjName, "\n")

prjDir <- "/home/users/skim/projects/collaboration/05_lineagetracing"

## sub-directories 
tabDir <- file.path(prjDir, "Tables"); tabDir
if(!dir.exists(tabDir)) dir.create(tabDir)

figDir <- file.path(prjDir, "Figures"); figDir
if(!dir.exists(figDir)) dir.create(figDir)

rdatDir <- file.path(prjDir, "Rdatasets"); rdatDir
if(!dir.exists(rdatDir)) dir.create(rdatDir)


## set working directory 
setwd( file.path(prjDir, "code"))



######################################################
## Sample information 
######################################################

dat <- readRDS(file.path (tabDir, "sinfo_LineageTracing_n267_sypark.RDS"))
dat$n_endomt_int <- round(dat$n_endomt)


## check count of cells per tissue 
tab.tissue <- sort(table(dat$tissue_id), decreasing=T)
table(tab.tissue)
# 1  2  3  4  5  7 
# 59 31 19 13  6  1 

dat0 <- dat; 
nrow(dat0) # 267 

##-----------------------------------------------------
## restrict data into tissues with at least two cells
dat <- dat0[dat0$tissue_id %in% names(tab.tissue)[tab.tissue>=2], ]
nrow(dat) # 208 


##----------------------------
## create tissueid_fac

tab.Source2 <- sort(table(dat$Source2), decreasing=T)
tab.Source_class <- sort(table(dat$Source_class), decreasing=T)

v.Source2 <- names(tab.Source2)
v.Source_class <- names(tab.Source_class)

dat$Source2_fac <- factor(dat$Source2, levels=v.Source2)
dat$Source_class_fac <- factor(dat$Source_class, levels=v.Source_class)

dat.ord <- dat[ order(dat$deadbody, 
                      dat$Source_class_fac, 
                      dat$Source2_fac, 
                      dat$tissue_id),]

v.tissue_id <- unique(dat.ord$tissue_id)

dat$tissue_id_fac <- factor(dat$tissue_id, levels=v.tissue_id)


##----------------------------------------------------
## remove age effect and plot 
##----------------------------------------------------

# fit.lmer <- lmer ( log2_n_endomt ~  age + (1 | tissue_id), data=dat)
# 
# summary(fit.lmer)
# anova(fit.lmer)
# coef(summary(fit.lmer))
# rand(fit.lmer)
# VarCorr(fit.lmer)
# dat$resid_noage <- residuals(fit.lmer)

fit.lm <- lm (log2_n_endomt ~ age, data=dat)
par(mfrow=c(2,2)); plot(fit.lm)
dat$resid_noage <- residuals(fit.lm)



##----------------------------------------------
## coloring

#library("scales")
#show_col(pal_npg("nrc")(10))
#show_col(pal_npg("nrc", alpha = 0.6)(10))

cols_deadbody <- pal_npg('nrc')(10)[c(1:5,9,7)]
cols_Source2 <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
#cols_Source2 <- distinctColorPalette(length(v.Source2))


## boxplot: baseline 
gp1 <- ggplot(dat, aes(x=tissue_id_fac, y=log2_n_endomt))+
  geom_boxplot(outlier.size=-1)+
  geom_point(aes(color=deadbody))+
  scale_color_manual(values = cols_deadbody)+
  ggtitle("log2 (n_endog_muts) ~ tissue_id") + 
  xlab("Skin tissues")+ylab("log2 (n_endog_muts)")  
  #theme_bw()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

gp1

## boxplot: resid ~ tissue id, color by deadbody 
gp2 <- ggplot(dat, aes(x=tissue_id_fac, y=resid_noage))+
  geom_boxplot(outlier.size=-1)+
  geom_point(aes(color=deadbody))+
  scale_color_manual(values = cols_deadbody)+
  ggtitle("Resid_noage ~ tissue_id") + 
  xlab("Skin tissues")+ylab("Residuals after regressing out Age")+
  geom_hline (yintercept=0, linetype="dashed", col="brown", size=1) + 
  geom_hline (yintercept=c(-0.5, 0.5), linetype="dashed", col=gray(0.5), size=1)  

gp2


## boxplot: resid ~ deadbody, color by deadbody
gp3 <- ggplot(dat, aes(x=deadbody, y=resid_noage))+
  geom_boxplot(outlier.size=-1)+
  geom_point(aes(color=deadbody))+ 
  geom_jitter(aes(colour = deadbody), width = 0.1) + 
  scale_color_manual(values = cols_deadbody) +
  ggtitle("Resid_noage ~ deadbody") + 
  xlab("Deadbody")+ylab("Residuals after regressing out Age")+
  geom_hline (yintercept=0, linetype="dashed", col="brown", size=1) + 
  geom_hline (yintercept=c(-0.5, 0.5), linetype="dashed", col=gray(0.5), size=1) 
gp3


## boxplot: resid ~ tissue id, color by source2
gp4 <- ggplot(dat, aes(x=tissue_id_fac, y=resid_noage))+
  geom_boxplot(outlier.size=-1)+
  geom_point(aes(color=Source2_fac))+
  scale_color_manual(values = cols_Source2)+
  ggtitle("Resid_noage ~ tissue_id") + 
  xlab("Skin tissues")+ylab("Residuals after regressing out Age")+
  geom_hline (yintercept=0, linetype="dashed", col="brown", size=1) + 
  geom_hline (yintercept=c(-0.5, 0.5), linetype="dashed", col=gray(0.5), size=1)  
gp4

## boxplot: resid ~ sosurce2, colour by source2
gp5 <- ggplot(dat, aes(x=Source2_fac, y=resid_noage))+
  geom_boxplot(outlier.size=-1)+
  geom_point(aes(color=Source2_fac))+ 
  geom_jitter(aes(colour = Source2_fac), width = 0.1) + 
  scale_color_manual(values = cols_Source2) +
  ggtitle("Resid_noage ~ Source2") + 
  xlab("Source2")+ylab("Residuals after regressing out Age")+
  geom_hline (yintercept=0, linetype="dashed", col="brown", size=1) + 
  geom_hline (yintercept=c(-0.5, 0.5), linetype="dashed", col=gray(0.5), size=1) 
gp5


png (file.path (figDir, "SUBSET_boxplot_check_donoreffect.png"), 
     width=1500, height=2000, res=130, pointsize=11)
## deadbody effect
grid.arrange(gp1, gp2, gp3, nrow=3)

dev.off()


png (file.path (figDir, "SUBSET_boxplot_check_Source2_tissue_effect.png"), 
     width=1500, height=1200, res=130, pointsize=11)
## other effect
grid.arrange(gp4, gp5, nrow=2)
dev.off()




##########################################################
## Try various modeling 
##########################################################

Mysummary.lm.detail <- function (fit) {
  res <- list(fit=fit,
              fit.summary=summary(fit),
              anova.test=anova(fit),
              coef.fixed=coef(summary(fit)))
  res
}


Mysummary.lmer.brief <- function (fit) {
  res <- list(fit.summary=summary(fit),
              space="---------------------------------------------------",
              rand.test=rand(fit))
  res
}


Mysummary.lmer.detail <- function (fit) {
  res <- list(fit=fit,
              fit.summary=summary(fit),
              anova.test=anova(fit),
              rand.test=rand(fit),
              coef.fixed=coef(summary(fit)),
              coef.rand=VarCorr(fit))
  res
}


Mysummary.glmer.detail <- function (fit) {
  res <- list(fit=fit,
              fit.summary=summary(fit),
              anova.test=anova(fit),
              coef.fixed=coef(summary(fit)),
              coef.rand=VarCorr(fit))
}



##==============================================
## linear models 
##==============================================

##----------------------------------------
## very basic: only age
fit0 <- lm ( log2_n_endomt ~ age, data=dat )
res.fit0 <- Mysummary.lm.detail (fit0)


##-------------------------------------------
## mixed effect modeling usig 'lmer'

## age + (1|deadbody)
fit1 <- lmer ( log2_n_endomt ~  age + (1|deadbody), data=dat)
res.fit1 <- Mysummary.lmer.detail (fit1) 

## age + (1|Source2)
fit2 <- lmer ( log2_n_endomt ~  age + (1|Source2_fac), data=dat)
res.fit2 <- Mysummary.lmer.detail (fit2) 

## age + (1|tissue_id)
fit3 <- lmer ( log2_n_endomt ~  age + (1|tissue_id_fac), data=dat)
res.fit3 <- Mysummary.lmer.detail (fit3) 

## age + (1|Source2) + (1|tissue_id)
fit4 <- lmer ( log2_n_endomt ~  age + (1|Source2_fac) + (1|tissue_id_fac:Source2_fac), data=dat)
res.fit4 <- Mysummary.lmer.detail (fit4) 



##==============================================
## generalized linear models 
##==============================================

##-------------------------------------------
## fixed effect model using 'glm'

## poisson 
gfit0 <- glm ( n_endomt_int ~  age , data=dat, family="poisson")
res.gfit0 <- Mysummary.lm.detail (gfit0)

## quasi-poisson allowing dispersion parameter
gfit0.quasi <- glm ( n_endomt_int ~  age , data=dat, family="quasipoisson")
res.gfit0.quasi <- Mysummary.lm.detail (gfit0.quasi)


##-------------------------------------------
## mixed effect modeling usig 'glmer'

if(0) {
  
  ## unstable..  & 'quasi-poisson' is not allowed for glmer
  gfit1 <- glmer ( n_endomt_int ~  age + (1 | deadbody), data=dat, family="poisson")
  # Warning message:
  #   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  #                  Model is nearly unidentifiable: very large eigenvalue
  #                - Rescale variables?
  
  gfit2 <- glmer ( n_endomt_int ~  age + (1 | tissue_id_fac), data=dat, family="poisson")
  # Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
  # Family: poisson  ( log )
  # Formula: n_endomt_int ~ age + (1 | tissue_id_fac)
  # Data: dat
  # AIC       BIC    logLik  deviance  df.resid 
  # 12584.201 12594.963 -6289.101 12578.201       264 
  # Random effects:
  #   Groups        Name        Std.Dev.
  # tissue_id_fac (Intercept) 0.2035  
  # Number of obs: 267, groups:  tissue_id_fac, 129
  # Fixed Effects:
  #   (Intercept)          age  
  # 6.54616      0.01381  
  # convergence code 0; 0 optimizer warnings; 2 lme4 warnings 
  
  
}



##-------------------------------------------
## try out variant component analysis 

#VCA::varPlot(form= resid_noage ~ deadbody, Data=dat)
#VCA::varPlot(form= resid_noage ~ Source2_fac, Data=dat)

png (file.path (figDir, "varPlot_VCA_resid-noage_Source2_tissue.png"), 
     width=1500, height=1000, res=130, pointsize=11)

VCA::varPlot(form= resid_noage ~ Source2_fac + tissue_id_fac, Data=dat)
abline(h=c(-0.5, 0, 0.5), col=gray(0.5), lwd=1.5, lty=c(2,1,2))

dev.off()



##===========================================
## fitting anova models 

vfit2 <- anovaVCA (resid_noage ~  Source2_fac, Data=dat)

vfit3 <- anovaVCA (resid_noage ~  tissue_id_fac, Data=dat)

vfit4 <- anovaVCA (resid_noage ~  Source2_fac + tissue_id_fac, Data=dat)


## portions explaind by each components
print(vfit2, 4)
print(vfit3, 4)
print(vfit4, 4)

#  > print(vfit2, 4)
# Result Variance Component Analysis:
#   -----------------------------------
#   
#   Name        DF      SS      MS     VC     %Total  SD     CV[%]               
# 1 total       228.313                0.1136 100     0.3371 -3593567926234798080
# 2 Source2_fac 24      5.6474  0.2353 0.0135 11.8686 0.1161 -1238014451154178048
# 3 error       242     24.2358 0.1001 0.1001 88.1314 0.3165 -3373581280953058816
# 
# Mean: 0 (N = 267) 
# 
# Experimental Design: unbalanced  |  Method: ANOVA
#
#-----------------------------------------------------------
#
# > print(vfit3, 4)
# Result Variance Component Analysis:
#   -----------------------------------
#   
#   Name          DF       SS      MS     VC     %Total SD     CV[%]               
# 1 total         243.6205                0.1126 100    0.3355 -3550295169312802816
# 2 tissue_id_fac 128      18.8114 0.147  0.0323 28.727 0.1798 -1902870557829203200
# 3 error         138      11.0718 0.0802 0.0802 71.273 0.2832 -2997278670626506752
# 
# Mean: 0 (N = 267) 
# 
# Experimental Design: unbalanced  |  Method: ANOVA
# 
#-----------------------------------------------------------
#
# > print(vfit4, 4)
# Result Variance Component Analysis:
#   -----------------------------------
#   
#   Name          DF       SS      MS     VC     %Total  SD     CV[%]               
# 1 total         220.0801                0.1135 100     0.3369 -3590485878239325696
# 2 Source2_fac   24       5.6474  0.2353 0.0099 8.7095  0.0994 -1059619704270729856
# 3 tissue_id_fac 105      13.5882 0.1294 0.0259 22.7965 0.1608 -1714301111007551488
# 4 error         137      10.6476 0.0777 0.0777 68.494  0.2788 -2971525975817055744
# 
# Mean: 0 (N = 267) 
# 
# Experimental Design: unbalanced  |  Method: ANOVA



##=========================================================
## confidence interval estimated from VCA? 
##=========================================================

#plotRandVar(vfit2, term="Source2_fac", mode="student")
inf.vfit2 <- VCAinference(vfit2,  VarVC=TRUE)
inf.vfit3 <- VCAinference(vfit3,  VarVC=TRUE)
inf.vfit4 <- VCAinference(vfit4,  VarVC=TRUE)


print(inf.vfit3)

######################################################
## Final summary 
######################################################


## significant age effect: 
res.fit0$fit.summary
res.fit1$fit.summary
res.gfit0.quasi$fit.summary ## reasonably similar t-value


## random effects
res.fit1$rand.test
cat("----------------------------------------------------\n")
res.fit2$rand.test
cat("----------------------------------------------------\n")
res.fit3$rand.test
cat("----------------------------------------------------\n")
res.fit4$rand.test

#   > res.fit1$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | deadbody)
# npar  logLik    AIC      LRT Df Pr(>Chisq)
# <none>            4 -74.785 157.57                       
# (1 | deadbody)    3 -74.801 155.60 0.030893  1     0.8605

# ----------------------------------------------------
#   > res.fit2$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | Source2_fac)
# npar  logLik    AIC    LRT Df Pr(>Chisq)    
# <none>               4 -63.864 135.73                         
# (1 | Source2_fac)    3 -74.801 155.60 21.875  1   2.91e-06 ***

# ----------------------------------------------------
#   > res.fit3$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | tissue_id_fac)
# npar  logLik    AIC    LRT Df Pr(>Chisq)    
# <none>                 4 -67.416 142.83                         
# (1 | tissue_id_fac)    3 -74.801 155.60 14.769  1  0.0001215 ***

# ----------------------------------------------------
#   > res.fit4$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | Source2_fac) + (1 | tissue_id_fac:Source2_fac)
# npar  logLik    AIC    LRT Df Pr(>Chisq)   
# <none>                             5 -62.476 134.95                        
# (1 | Source2_fac)                  4 -66.955 141.91 8.9585  1   0.002762 **
# (1 | tissue_id_fac:Source2_fac)    4 -63.864 135.73 2.7760  1   0.095684 . 



