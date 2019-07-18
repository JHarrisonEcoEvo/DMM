rm(list=ls())
options(scipen = 99)

library(DESeq2)
library(edgeR)
library(rstan)
library(exactRankTests)
library(nlme)
library(ggplot2)
library(ALDEx2)
set.seed(124)

##################
# Begin analyses #
##################

#wrangle
dat <- read.table("./data/rosen_mincount10_maxee2_trim200_results_forpaper/rosen_mincount10_maxee2_trim200.otu_table.99.denovo.rdp_assigned.paper_samples.txt",
                  fill = T)
metadat <- read.csv("./data/patient_clinical_metadata.csv")

#The two groups being compared are delineated in metadat$mbs_consolidated (aspirators vs. normal)

#Remove the subjects that were not measured
metadat <- metadat[which(!(metadat$mbs_consolidated == "")),]

#split off the designator for sample location from the otu table
#bal = bronchoalveolar lavage
#g = gastric fluid
#T = swab
#TI = ?

#Remove subjects ending in FI or TI
dat <- dat[-grep("TI$", row.names(dat)),]
dat <- dat[-grep("F1$", row.names(dat)),]
dat <- dat[-grep("F1T$", row.names(dat)),]
dat <- dat[-grep("SI$", row.names(dat)),]

location <- gsub("\\d+-\\d+-\\d+([GTISB]*)$", "\\1",
                 row.names(dat))

location <- gsub("\\d+-\\d+-\\d+F1([GTISB]*)$", "\\1",
                 location)   

dat$location <- location

#Now remove location from the sample names and merge with metadata
samps <- gsub("(\\d+-\\d+-\\d+)[GTISB]*$", "\\1",
              row.names(dat))
dat$samps <- samps

dat <- merge(metadat, dat, by.x = "subject_id", by.y="samps")

dim(dat)
table(dat$location, dat$mbs_consolidated,  dat$gender_all)
#might be worth doing for either sex separately

#order otu table by metadata
#first subset to just B samples, will need to run for each
baldat <- dat[dat$location == "B",]
baldat <- baldat[order(baldat$mbs_consolidated),]

#get grouping indices
getgroups <- function(x){
  starts = c(min(which(x == unique(x)[1])),
             min(which(x == unique(x)[2]))) 
  ends = c(max(which(x == unique(x)[1])),
           max(which(x == unique(x)[2])))
  
  return(list(
    starts=starts,
    ends = ends))
}
getgroups(baldat$mbs_consolidated)

#compile model
DM <- rstan::stan_model("./DM.stan", 
                        model_name="DM")

#run model
fitstan_HMC <-rstan::sampling(DM, 
                             data=list("datamatrix"=baldat[,34:4039]+1, 
                                       "nreps"=nrow(baldat[,34:4039]), 
                                       "notus"=ncol(baldat[,34:4039]),
                                       "N"=2, 
                                       "start" = getgroups(baldat$mbs_consolidated)$starts,
                                       "end" = getgroups(baldat$mbs_consolidated)$ends), 
                             chains=2, 
                             control = list(max_treedepth = 15),
                             warmup=1000, 
                             iter=2000, 
                             thin=2, 
                             algorithm="NUTS", 
                             cores=2, 
                             seed=123,
                             pars<-c("pi")
)

save.image(file = "./reanalysisLung.RData")

# Perform analysis with VB, if desired.
# fitstan_VB <-rstan::vb(DM,
#                              data=list("datamatrix"=baldat[,34:4039]+1,
#                                        "nreps"=nrow(baldat[,34:4039]),
#                                        "notus"=ncol(baldat[,34:4039]),
#                                        "N"=2,
#                                        "start" = getgroups(baldat$mbs_consolidated)$starts,
#                                        "end" = getgroups(baldat$mbs_consolidated)$ends),
#                             # chains=2,
#                              #control = list(max_treedepth = 15),
#                             # warmup=1000,
#                             # iter=2000,
#                             # thin=2,
#                              algorithm="meanfield",
#                              #cores=2,
#                              seed=123,
#                              pars<-c("pi")
# )
# 
# 
# est.pi<-extract(fitstan_VB,"pi")
# 
# #Differential abundance testing.
# 
# diffs_VB <- est.pi$pi[,1,] - est.pi$pi[,2,]
# 
# calc_certain_diffs <- function(mcmc_of_diffs){
#   positives <- vector()
#   negatives <- vector()
#   
#   for(i in 1:dim(mcmc_of_diffs)[2]){
#     perc <- length(which(mcmc_of_diffs[,i] > 0 ))/ length(mcmc_of_diffs[,i])
#     if(perc >= 0.95 | perc <= 0.05){
#       positives <- c(positives, i)
#     }else{
#       negatives <- c(negatives, i)
#     }
#   }
#   return(list(positives = positives,
#               negatives = negatives))
# }
# 
# outVB <- calc_certain_diffs(diffs_VB)
# 
# 
# #Bring in HMC data and compare to VB. 
# #Note that HMC implementation takes time time to run.
# load(file = "./reanalysisLung.RData")
# # shinystan::launch_shinystan(fitstan_VB, rstudio=T) 
# 
# est.pi<-extract(fitstan_VB,"pi")
# 
# #do differential abundance tests.
# 
# diffs_VB <- est.pi$pi[,1,] - est.pi$pi[,2,]
# 
# calc_certain_diffs <- function(mcmc_of_diffs){
#   positives <- vector()
#   negatives <- vector()
#   
#   for(i in 1:dim(mcmc_of_diffs)[2]){
#     perc <- length(which(mcmc_of_diffs[,i] > 0 ))/ length(mcmc_of_diffs[,i])
#     if(perc >= 0.95 | perc <= 0.05){
#       positives <- c(positives, i)
#     }else{
#       negatives <- c(negatives, i)
#     }
#   }
#   return(list(positives = positives,
#               negatives = negatives))
# }
# 
# outNUTS <- calc_certain_diffs(diffs_VB)
# length(outNUTS$positives)
# 
# setdiff(outVB$positives, outNUTS$positives)
# setdiff(outNUTS$positives, outVB$positives)
# length(outVB$positives)
# length(outNUTS$positives)
# 
# rev(sort(apply(diffs_VB, 2, mean)))
# sort(apply(diffs_VB, 2, mean))
# 
# which(apply(diffs_VB, 2, mean) > 0.01)
# names(baldat[,34:4039])[1684]
# 
# 
# 
# #Make a figure showing what changed. 
# forfig <- diffs_VB[,outNUTS$positives]
# 
# resultsPlotter <- function(x) {
#   plot(
#     NULL,
#     ylim = c(min(apply(
#       x, 2, FUN = quantile, probs = c(0.025)
#     )),
#     max(apply(
#       x, 2, FUN = quantile, probs = c(0.975)
#     ))),
#     xlim = c(0, dim(x)[2]),
#     xaxt = "n",
#     yaxt = "n",
#     frame = F,
#     ylab = "",
#     xlab = ""
#   )
#   axis(
#     side = 1,
#     at = c(seq(0, dim(x)[2], by = 10), dim(x)[2]),
#     labels = c(seq(0, dim(x)[2], by = 10), dim(x)[2])
#   )
#   
#   axis(
#     side = 2,
#     at = c(seq(min(apply(x, 2, FUN = quantile, probs = c(0.025))), 
#                max(apply(x, 2, FUN = quantile, probs = c(0.975))), 
#                by = round(max(apply(x, 2, FUN = quantile, probs = c(0.975)))/5, digits = 3))),
#     labels = round(abs(c(seq(min(apply(x, 2, FUN = quantile, probs = c(0.025))),  
#                              max(apply(x, 2, FUN = quantile, probs = c(0.975))), 
#                              by = round(max(apply(x, 2, FUN = quantile, probs = c(0.975)))/5, digits = 3)))), digits=3),
#     las = 2
#   )  
#   points(apply(x, 2, FUN = mean),
#          pch = 18)
#   
#   segments(
#     y0 = apply(x, 2, FUN = quantile, probs = c(0.025)),
#     x0 = seq(1, dim(x)[2]),
#     y1 = apply(x, 2, FUN = quantile, probs = c(0.975)),
#     x1 = seq(1, dim(x)[2])
#   )
#   abline(h = 0,
#          lty = 3)
# }
# pdf(file = "./visuals/lungHMC.pdf", width = 7, height=7)
# par(oma = c(0,2,0,0))
# resultsPlotter(forfig[,order(apply(forfig, 2, FUN = mean))])
# text(x = 40, 
#      y = 0.018,
#      "aspirators")
# text(x = 40, 
#      y = -0.005,
#      "non-aspirators")
# mtext(side = 2, 
#       "Effect size (proportion total reads)",
#       line = 4)
# mtext(side = 1, 
#       "Taxon",
#       line = 2)
# dev.off()