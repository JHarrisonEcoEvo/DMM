#J. Harrison
library(xtable)
rm(list=ls())
set.seed(666)
options(scipen = 99)

dat <- read.table("./data/allv4.csv", 
                  stringsAsFactors = F)

#remove the unwanted header names, after improving the names of the data.frame
names(dat) <- dat[3,]

dat <- dat[-grep("notus", dat$notus),]

#SANITY CHECK
dim(dat) #Check that dimensions are correct. 

dat$combo <- paste(dat$modelType, 
                   dat$notus, 
                   dat$nreads,
                   dat$nreps,
                   dat$precision_theta,
                   dat$num_taxa_different,
                   dat$effectsize,
                   #dat$simRep,
                   sep = "")

table(dat$notus)
table(dat$nreads)
table(dat$modelType)

# 
# #SANITY CHECK
# #check that everything ran. 
# params <- read.csv("serialRunParamsv2.csv")
# sanityCombo <- paste(params[,1], 
#                      params[,2], 
#                      params[,3], 
#                      params[,4], 
#                      params[,5], 
#                      params[,6], 
#                      params[,7]
#                      #params[,8]
#                      , sep = "")
# 
# #Note the !
# notdone <- paste(params[,1], 
#                  params[,2], 
#                  params[,3], 
#                  params[,4], 
#                  params[,5], 
#                  params[,6], 
#                  params[,7]
#                  #params[,8]
#                  ,sep = ":")[which(!(sanityCombo %in% dat$combo))]
# 
# redo <- data.frame(matrix(unlist(strsplit(notdone, ":")), 
#                           byrow = T,
#                           ncol = 8))
# dim(redo)
# write.csv(redo, file = "redo.csv", row.names = F)

#REMAKE the combo variable without simrep, so that 
#averaging across replicates works properly
dat <- dat[,-which(names(dat) == "combo")]
dim(dat)
combo <- paste(dat$modelType, 
               dat$notus, 
               dat$nreads,
               dat$nreps,
               dat$precision_theta,
               dat$num_taxa_different,
               dat$effectsize,
               sep = "")
dat <- data.frame(combo,dat)
names(dat)[1] <- "combo"

#Convert everything to numbers, to clean up code below. 

for(i in 3:length(dat)){
  dat[,i] <- as.numeric(dat[,i])
}

#Convert true positives to proportions. 

truepPropsAll <- dat[, grep("truep_\\w+_all", names(dat))] / dat$numDiffering
truepPropsabund <- dat[, grep("truep_\\w+_abund", names(dat))] / dat$numDiffering_abund
truepPropsmed <- dat[, grep("truep_\\w+_med$", names(dat))] / dat$numDiffering_med
truepPropsrare <- dat[, grep("truep_\\w+_rare", names(dat))] / dat$numDiffering_rare
truepPropsmedhi <- dat[, grep("truep_\\w+_medhi", names(dat))] / dat$numDiffering_rare

#Calculate false positive rates
dat$fpr_diff_all <- dat$falsep_diff_all/(dat$falsep_diff_all + dat$truen_diff_all) 
dat$fpr_wilcox_all <- dat$falsep_wilcoxBH_FDR_all/(dat$falsep_wilcoxBH_FDR_all + dat$truen_wilcoxBH_FDR_all) 
dat$fpr_deseq_all <- dat$falsep_deseq_all/(dat$falsep_deseq_all + dat$truen_deseq_all) 
dat$fpr_edger_all <- dat$falsep_edgeR_all/(dat$falsep_edgeR_all + dat$truen_edgeR_all) 
dat$fpr_ancom_all <- dat$falsep_ancom_all/(dat$falsep_ancom_all + dat$truen_ancom_all) 
dat$fpr_stan_all <- dat$falsep_stanHMC_all/(dat$falsep_stanHMC_all + dat$truen_stanHMC_all) 
dat$fpr_stanvb_all <- dat$falsep_stanVB_all/(dat$falsep_stanVB_all + dat$truen_stanVB_all) 
dat$fpr_aldexT_all <- dat$falsep_aldexT_all/(dat$falsep_aldexT_all + dat$truen_aldexT_all) 
dat$fpr_jags_clr_all <- dat$falsep_jags_clr_all/(dat$falsep_jags_clr_all + dat$truen_jags_clr_all) 

dat <- data.frame(dat, truepPropsAll, truepPropsabund, truepPropsmed, truepPropsrare, truepPropsmedhi)
out <- dat
####################
# Perform analysis #
####################
#Make Uniform the reference level
out$modelType <- as.factor(out$modelType)
out$modelType <- relevel(out$modelType, ref = "uniform")

regplotr <- function(reg){
  par(mfrow = c(2,2))
  plot(reg, which = c(1,2,3,4))
}

results <- data.frame(
                  matrix(ncol = 9,
                  nrow = 5))

k <- 1
for(i in grep("all.1$", names(out))){
  reg <- lm(out[,i] ~ out$modelType
    + out$notus
    + out$nreads
    + out$nreps
    + out$precision_theta
    )
  
  #Calculate CIs in beta coefficients. 
  his <- round(summary(reg)$coefficients[,1] + 1.96*summary(reg)$coefficients[,2], 
               2)
  lows <- round(summary(reg)$coefficients[,1] - 1.96*summary(reg)$coefficients[,2],
              2)
  
  stars <- vector()
  stars[which(round(summary(reg)$coefficients[,4], 2) <= 0.1)] <- "*"
  stars[which(round(summary(reg)$coefficients[,4], 2) <= 0.05)] <- "**"
  stars[which(round(summary(reg)$coefficients[,4], 2) <= 0.01)] <- "***"
  stars[which(round(summary(reg)$coefficients[,4], 2) > 0.1)] <- ""
  
  
  #output results
  results[k,1] <- names(out)[i]
    
  results[k,2:8] <- paste(round(summary(reg)$coefficients[,1], 2),
    paste("(",
        paste(his, lows, sep = ", "),
        ")",
        stars,
        sep = "")
  )

  regplotr(reg)
  k <- k + 1
}

names(results) <- c("Method",
                    "Intercept",
                    "Pareto (shape = 0.7)",
                    "Pareto (shape = 4)",
                    "Num. of features",
                    "Num. of observations",
                    "Num. of replicates",
                    "Noise")

#Rewrite so it looks nice in latex
results[,1] <- c("DMM Gibbs",
                    "Wilcoxon",
                    "DESeq2",
                    "edgeR",
                    "ANCOM",
                 "DMM HMC",
                 "DMM VB",
                 "DMM Gibbs CLR",
                 "Aldex T-test")

print(xtable(results, type = "latex"),
      file = "lm_results.tex",
      floating.environment='sidewaystable',
      include.rownames=FALSE)

#sanity check

lm(out$truep_ancom_all.1 ~ out$modelType
   + out$notus
   + out$nreads
   + out$nreps
   + out$precision_theta
)

lm(out$truep_diff_all.1 ~ out$modelType
   + out$notus
   + out$nreads
   + out$nreps
   + out$precision_theta
)

#fpr

results <- data.frame(
  matrix(ncol = 9,
         nrow = 5))

k <- 1
for(i in grep("fpr[A-Za-z_]*all", names(out))){
  reg <- lm(out[,i] ~ out$modelType
            + out$notus
            + out$nreads
            + out$nreps
            + out$precision_theta
  )
  
  #Calculate CIs in beta coefficients. 
  his <- round(summary(reg)$coefficients[,1] + 1.96*summary(reg)$coefficients[,2], 
               2)
  lows <- round(summary(reg)$coefficients[,1] - 1.96*summary(reg)$coefficients[,2],
                2)
  
  stars <- vector()
  stars[which(round(summary(reg)$coefficients[,4], 2) <= 0.1)] <- "*"
  stars[which(round(summary(reg)$coefficients[,4], 2) <= 0.05)] <- "**"
  stars[which(round(summary(reg)$coefficients[,4], 2) <= 0.01)] <- "***"
  stars[which(round(summary(reg)$coefficients[,4], 2) > 0.1)] <- ""
  
  
  #output results
  results[k,1] <- names(out)[i]
  
  results[k,2:8] <- paste(round(summary(reg)$coefficients[,1], 2),
                          paste("(",
                                paste(his, lows, sep = ", "),
                                ")",
                                stars,
                                sep = "")
  )
  
  #regplotr(reg)
  k <- k + 1
}

names(results) <- c("Method",
                    "Intercept",
                    "Pareto (shape = 0.7)",
                    "Pareto (shape = 4)",
                    "Num. of features",
                    "Num. of observations",
                    "Num. of replicates",
                    "Noise")

#Rewrite so it looks nice in latex
results[,1] <- c("DMM Gibbs",
                 "Wilcoxon",
                 "DESeq2",
                 "edgeR",
                 "ANCOM",
                 "DMM HMC",
                 "DMM VB",
                 "DMM Gibbs CLR",
                 "Aldex T-test")

print(xtable(results, type = "latex"),
      file = "lm_results_fdr.tex",
      floating.environment='table',
      include.rownames=FALSE)

#sanity check
re<-lm(out$fpr_ancom_all ~ out$modelType
   + out$notus
   + out$nreads
   + out$nreps
   + out$precision_theta
);summary(re)
re<- lm(out$fpr_diff_all ~ out$modelType
   + out$notus
   + out$nreads
   + out$nreps
   + out$precision_theta
);summary(re)
