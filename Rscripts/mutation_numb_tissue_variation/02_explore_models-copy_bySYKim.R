##
## Model mutation counts 
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

nrow(dat) ## 267 
 
## main IDs 
dat$deadbody  ## donor
dat$age  ## age
dat$gender  ## gender
dat$batch ## sequencing batch

dat$sample_id # cell id

dat$tissue_id # tissue id
dat$Source2  # anatomy id1
dat$Source_class #anatomy id2

dat$n_endomt  ## endogeneous mutation counts
dat$log2_n_endomt ## log2 scale



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
## endogeneous mutations show linear relationship with deadbody 
##----------------------------------------------------

png (file.path (figDir, "Assoc_endoMuts_Age_donor.png"), 
     width=1000, height=1000, res=130, pointsize=11)

par(mfrow=c(2,2))
# original scale
plot(dat$age, dat$n_endomt)
abline(lm(dat$n_endomt~dat$age), col="red")
b <- boxplot(dat$n_endomt ~ dat$deadbody)
points(jitter(as.numeric(dat$deadbody), 1), dat$n_endomt, col=gray(0.5), cex=0.5)

## log2 scale
plot(dat$age, dat$log2_n_endomt)
abline(lm(dat$log2_n_endomt~dat$age), col="red")
b <- boxplot(dat$log2_n_endomt ~ dat$deadbody)
points(jitter(as.numeric(dat$deadbody), 1), dat$log2_n_endomt, col=gray(0.5), cex=0.5)

dev.off()



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



png (file.path (figDir, "boxplot_check_donoreffect.png"), 
     width=1500, height=2000, res=130, pointsize=11)
## deadbody effect
grid.arrange(gp1, gp2, gp3, nrow=3)

dev.off()


png (file.path (figDir, "boxplot_check_Source2_tissue_effect.png"), 
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

## age + (1|Source2) + (1|tissue_id:Source2_fac)
fit4 <- lmer ( log2_n_endomt ~  age + (1|Source2_fac) + (1|tissue_id_fac:Source2_fac), data=dat)
res.fit4 <- Mysummary.lmer.detail (fit4) 

## age + (1|Source2) + (1|tissue_id)
fit5 <- lmer ( log2_n_endomt ~  age + (1|Source2_fac) + (1|tissue_id_fac), data=dat)
res.fit5 <- Mysummary.lmer.detail (fit5) 

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
summary(fit4)
res.fit5$rand.test
summary(fit5)

# > res.fit1$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | deadbody)
# npar  logLik    AIC     LRT Df Pr(>Chisq)
# <none>            4 -94.864 197.73                      
# (1 | deadbody)    3 -95.102 196.20 0.47557  1     0.4904

# ----------------------------------------------------
#   > res.fit2$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | Source2_fac)
# npar  logLik   AIC    LRT Df Pr(>Chisq)    
# <none>               4 -86.800 181.6                         
# (1 | Source2_fac)    3 -95.102 196.2 16.605  1  4.603e-05 ***

# ----------------------------------------------------
#   > res.fit3$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | tissue_id_fac)
# npar  logLik    AIC    LRT Df Pr(>Chisq)    
# <none>                 4 -87.970 183.94                         
# (1 | tissue_id_fac)    3 -95.102 196.20 14.264  1  0.0001589 ***

# ----------------------------------------------------
#   > res.fit4$rand.test
# ANOVA-like table for random-effects: Single term deletions
# 
# Model:
#   log2_n_endomt ~ age + (1 | Source2_fac) + (1 | tissue_id_fac:Source2_fac)
# npar  logLik    AIC    LRT Df Pr(>Chisq)  
# <none>                             5 -84.290 178.58                       
# (1 | Source2_fac)                  4 -87.486 182.97 6.3922  1    0.01146 *
# (1 | tissue_id_fac:Source2_fac)    4 -86.800 181.60 5.0183  1    0.02508 *






