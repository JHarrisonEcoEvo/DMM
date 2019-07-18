#J. Harrison
rm(list=ls())
set.seed(666)
options(scipen = 99)

dat <- read.table("./data/allv4.csv", 
                  stringsAsFactors = F)
#remove the unwanted header names, after improving the names of the data.frame
names(dat) <- dat[1,]

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


#SANITY CHECK
#check that everything ran. 
params <- read.csv("serialRunParamsv2.csv")
sanityCombo <- paste(params[,1], 
                     params[,2], 
                     params[,3], 
                     params[,4], 
                     params[,5], 
                     params[,6], 
                     params[,7]
                     #params[,8]
                     , sep = "")

#Note the !
notdone <- paste(params[,1], 
                 params[,2], 
                 params[,3], 
                 params[,4], 
                 params[,5], 
                 params[,6], 
                 params[,7]
                 #params[,8]
                ,sep = ":")[which(!(sanityCombo %in% dat$combo))]

redo <- data.frame(matrix(unlist(strsplit(notdone, ":")), 
                          byrow = T,
                          ncol = 7))
dim(redo)
write.csv(redo, file = "redo.csv", row.names = F)

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

#Combine model performance among reps by averaging, do this for all model types not just the DMM
# out <- data.frame(matrix(ncol = 293, 
#                          nrow = length(unique(dat$combo))))
# k <- 1
# for(i in unique(dat$combo)){
#   d <- subset(dat, combo == i)
#   out[k,] <- c(d[1, c(1:9)], apply(d[,10:187], 2, FUN = mean))
#   k <- k + 1
# }
# 
# names(out) <- names(dat)[c(1:9, 10:187)]

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

colz <- c("coral2",#RColorBrewer::brewer.pal(n = 5, name = "Accent")[1], 
          "gold3", 
          "cyan4") #c("darkorange","darkorchid","cyan3")
########
# Plot #
########

#what are the outlier rmse for vb?
dat[rev(order(dat$RMSE))[1:15],1:7]

pdf(file = "./visuals/rmse.pdf",
    width = 8,
    height = 8)
#RMSE 
par(mfrow = c(2,2),
    mar = c(2,3,1,1),
    oma = c(6,5,1,1))
stripchart(dat$RMSE ~ 
             dat$modelType,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = #c(add.alpha(colz, alpha = 0.9),
             add.alpha(colz, alpha = 0.5),
           # add.alpha(colz, alpha = 0.2)),
           cex=1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           ylim=c(0,0.01),
           xlim = c(0,3.5),
           frame.plot = FALSE,
           xpd = NA
)
axis(side = 2, 
     at = c(seq(0,0.01, by =0.005)),
     labels = c(seq(0,0.01, by = 0.005)),
     las = 2, 
     cex.axis = 1.8)

boxplot(dat$RMSE~ 
          dat$modelType,
        las=2,
        yaxt="n",
        xaxt="n",
        outline=F,
        add=T,
        col =# c(add.alpha(colz, alpha = 0.9),
          add.alpha(colz, alpha = 0.5),
        # add.alpha(colz, alpha = 0.2)),
        xpd = NA,
        xlim = c(0,3.5),
        frame.plot = FALSE
)

text("(a)",
     xpd = NA, 
     y = 0.01,
     x = 0.5,
     cex = 2)
#Plot b
mtext("RMSE", line= 4.7, side = 2, cex = 2, adj = -.14)
stripchart(dat$rmseAbund ~ 
             dat$modelType,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = 
             add.alpha(colz, alpha = 0.5)[1:2],
           cex=1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           ylim=c(0,0.01),
           xlim = c(0,2.5),
           frame.plot = FALSE,
           xpd = NA
)
axis(side = 2, 
     at = c(seq(0,0.01, by =0.005)),
     labels = c(seq(0,0.01, by = 0.005)),
     las = 2, 
     cex.axis = 1.8)

boxplot(dat$rmseAbund ~ 
          dat$modelType,
        las=2,
        yaxt="n",
        xaxt="n",
        outline=F,
        add=T,
        col = add.alpha(colz, alpha = 0.5)[1:2],
        xpd = NA,
        xlim = c(0,2.5),
        frame.plot = FALSE
)
text("(b)",
     xpd = NA, 
     y = 0.01,
     x = 0.4,
     cex = 2)
#Plot c

stripchart(dat$rmseMed ~ 
             dat$modelType,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = 
             add.alpha(colz, alpha = 0.5)[1:3],
           cex=1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           ylim=c(0,0.01),
           xlim = c(0,3.5),
           frame.plot = FALSE,
           xpd = NA
)
axis(side = 2, 
     at = c(seq(0,0.01, by =0.005)),
     labels = c(seq(0,0.01, by = 0.005)),
     las = 2, 
     cex.axis = 1.8)


boxplot(dat$rmseMed ~ 
          dat$modelType,
        las=2,
        yaxt="n",
        xaxt="n",
        outline=F,
        add=T,
        col = add.alpha(colz, alpha = 0.5)[1:3],
        xpd = NA,
        xlim = c(0,3.5),
        frame.plot = FALSE
)
text("(c)",
     xpd = NA, 
     y = 0.01,
     x = 0.5,
     cex = 2)
#Plot d

stripchart(dat$rmseRare ~ 
             dat$modelType,
           vertical = TRUE, 
           data = dat, 
           method = "jitter", 
           pch = 20, 
           col = 
             add.alpha(colz, alpha = 0.5)[1:2],
           cex=1.5,
           xaxt="n",
           yaxt="n",
           ylab="",
           ylim=c(0,0.01),
           xlim = c(0,2.5),
           frame.plot = FALSE,
           xpd = NA
)
axis(side = 2, 
     at = c(seq(0,0.01, by =0.005)),
     labels = c(seq(0,0.01, by = 0.005)),
     las = 2, 
     cex.axis = 1.8)


boxplot(dat$rmseRare ~ 
          dat$modelType,
        las=2,
        yaxt="n",
        xaxt="n",
        outline=F,
        add=T,
        col = add.alpha(colz, alpha = 0.5)[1:2],
        xpd = NA,
        xlim = c(0,2.5),
        frame.plot = FALSE
)
text("(d)",
     xpd = NA, 
     y = 0.01,
     x = 0.4,
     cex = 2)

#Legend
rect(xleft = -2.95,
     xright = -2.45,
     ytop = -0.0010,
     ybottom = -0.002,
     border = NA,
     col = add.alpha(colz[1], alpha = 0.9),
     xpd = NA)

text("Pareto (0.7)",
     x = -2.5,
     y = -0.0015,
     xpd = NA,
     pos = 4,
     cex = 1.5)

rect(xleft = -1.15+0.25,
     xright = -0.65+0.25,
     ytop = -0.0010,
     ybottom = -0.002,
     border = NA,
     col = add.alpha(colz[2], alpha = 0.9),
     xpd = NA)

text("Pareto (4)",
     x = -0.7+0.25,
     y = -0.0015,
     xpd = NA,
     pos = 4,
     cex = 1.5)

rect(xleft = 0.65+0.25,
     xright = 1.15+0.25,
     ytop = -0.0010,
     ybottom = -0.002,
     border = NA,
     col = add.alpha(colz[3], alpha = 0.9),
     xpd = NA)

text("Equal",
     x = 1.1+0.25,
     y = -0.0015,
     xpd = NA,
     pos = 4,
     cex = 1.5)
dev.off()


##################################
# RMSE of three model types #
##################################
pdf(width=8,height = 8, file = "./visuals/rmse_threeMCMC.pdf")
par(oma=c(3,4,2,2), mfrow =c(1,1))

stripchart(list(dat$RMSE, 
                dat$rmse_hmc,
                dat$rmse_vb),
           vertical = TRUE, 
           method = "jitter", 
           pch = 20, 
           col = 
             add.alpha(c("coral3","darkseagreen4","darkslateblue"), alpha = 0.5),
           cex=1.5,
           xaxt="n",
           yaxt="n",
           #group.names=c("Gibbs", "HMC", "VB"),
           #cex.lab = 2,
           ylim=c(0,0.003),
           xlim = c(0.5,3.5),
           frame.plot = FALSE,
           xpd = NA
)
mtext(side=2, "RMSE", line =5.5, cex = 2.5)
axis(side = 2, 
     at = c(seq(0,0.003, by = 0.001)),
     labels = c(seq(0,0.003, by = 0.001)),
     las = 2, 
     cex.axis = 1.8)
text( c("MCMC", "HMC", "VI"),
      x = c(1,2,3),
      y = -0.0005, 
     cex = 1.9, 
     xpd = NA)

boxplot(list(dat$RMSE, 
             dat$rmse_hmc,
             dat$rmse_vb),
        las=2,
        yaxt="n",
        xaxt="n",
        outline=F,
        add=T,
        col =
        add.alpha(c("coral3","darkseagreen4","darkslateblue"), alpha = 0.5),
        xpd = NA,
        frame.plot = FALSE
)

dev.off()

##################################
# time of three model types #
##################################
pdf(width=8,height = 8, file = "./visuals/time_threeMCMC.pdf")
par(oma=c(3,4,2,2), mfrow =c(1,1))

stripchart(list(dat$dmm_time,
                dat$stanHMCtime,
                dat$stanVBtime),
           vertical = TRUE,
           method = "jitter",
           pch = 20,
           col =
             add.alpha(colz, alpha = 0.5)[1:2],
           cex=1.5,
           xaxt="n",
           #yaxt="n",
           #group.names=c("Gibbs", "HMC", "VB"),
           #cex.lab = 2,
          # ylim=c(0,0.04),
           xlim = c(0.5,3.5),
           frame.plot = FALSE,
           xpd = NA
)
mtext(side=2, "Time (s)", line =5, cex = 2)
# axis(side = 2,
#      at = c(seq(0,0.04, by = 0.01)),
#      labels = c(seq(0,0.04, by = 0.01)),
#      las = 2,
#      cex.axis = 1.8)
text( c("MCMC", "HMC", "VB"),
      x = c(1,2,3),
      y = -0.005,
      cex = 1.5,
      xpd = NA)

boxplot(list(dat$dmm_time,
             dat$stanHMCtime,
             dat$stanVBtime),
        las=2,
        yaxt="n",
        xaxt="n",
        outline=F,
        add=T,
        col = add.alpha(colz, alpha = 0.5)[1:2],
        xpd = NA,
        frame.plot = FALSE
)

dev.off()

##################################
# hdi inclusivity  #
##################################

reg <- lm(dat$perc_overlap~dat$actualOTUs + dat$modelType + dat$nreps + dat$nreads + dat$precision_theta)
summary(reg)

boxplot(dat$perc_overlap~dat$modelType)
boxplot(dat$RMSE~dat$modelType)
boxplot(dat$truth_in_hdi_hmc~dat$modelType)
boxplot(dat$truep_stanVB_all.1~dat$modelType)
boxplot(dat$truep_diff_all.1~dat$modelType)
boxplot(dat$truep_stanHMC_all.1~dat$modelType)

pdf(width=8,height = 8, file = "./visuals/hdi_threeMCMC.pdf")
par(oma=c(3,4,2,2), mfrow =c(1,1))


stripchart(list(dat$perc_overlap, 
                dat$truth_in_hdi_hmc,
                dat$truth_in_hdi_vb),
           vertical = TRUE, 
           method = "jitter", 
           pch = pointtype, 
           col = 
             add.alpha(colz, alpha = 0.5)[1:2],
           cex=1.5,
           xaxt="n",
           #yaxt="n",
           #group.names=c("Gibbs", "HMC", "VB"),
           #cex.lab = 2,
            ylim=c(0,1),
           xlim = c(0.5,3.5),
           frame.plot = FALSE,
           xpd = NA
)


mtext(side=2, "Proportion of truth within HDI", line =5, cex = 2)
# axis(side = 2, 
#      at = c(seq(0,0.04, by = 0.01)),
#      labels = c(seq(0,0.04, by = 0.01)),
#      las = 2, 
#      cex.axis = 1.8)
text( c("MCMC", "HMC", "VB"),
      x = c(1,2,3),
      y = -0.04, 
      cex = 1.5, 
      xpd = NA)

boxplot(list(dat$perc_overlap, 
             dat$truth_in_hdi_hmc,
             dat$truth_in_hdi_vb),
        las=2,
        yaxt="n",
        xaxt="n",
        outline=F,
        add=T,
        col = add.alpha(colz, alpha = 0.5)[1:2],
        xpd = NA,
        frame.plot = FALSE
)

dev.off()


##################################
# truep of three model types by rank abundance #
##################################
pdf(width=8,height = 8, file = "./visuals/perfhmcvbgibbs_rankabundance.pdf")
par(oma=c(3,3,0,0),mar=c(3,3,3,3), mfrow =c(3,2))

boxplot(dat$perc_overlap~dat$modelType,
        main = "", 
        las = 2,
        names = c("","",""), 
        col=add.alpha(colz, alpha = 0.5),
        outline=F,
        frame = F)
mtext(side=3, "Proportion of features with\n truth in HDI")
mtext(side=2, "MCMC",
      line=3)

boxplot(dat$truep_diff_all.1~dat$modelType,
        main = "",
        las = 2,
        names = c("","",""), 
        col=add.alpha(colz, alpha = 0.5),
        outline=F,
        frame = F)
mtext(side=3, "Proportion of true positives recovered")

boxplot(dat$truth_in_hdi_hmc~dat$modelType,
        main = "", 
        las = 2,
        names = c("","",""), 
        col=add.alpha(colz, alpha = 0.5),
        outline=F,
        frame = F)
mtext(side=2, "HMC",
      line=3)

boxplot(dat$truep_stanHMC_all.1~dat$modelType,
        main = "", 
        las = 2,
        names = c("","",""), 
        col=add.alpha(colz, alpha = 0.5),
        outline=F,
        frame = F)

boxplot(dat$truth_in_hdi_vb~dat$modelType,
        main = "",
        las = 2, 
        col=add.alpha(colz, alpha = 0.5),
        outline=F,
        names=c("Pareto, 0.7", "Pareto, 4", "Equal"),
        frame = F)
mtext(side=2, "VB",
      line=3)

boxplot(dat$truep_stanVB_all.1~dat$modelType,
        main = "",
        las = 2, 
        col=add.alpha(colz, alpha = 0.5),
        outline=F,
        names=c("Pareto, 0.7", "Pareto, 4", "Equal"),
        frame = F)

dev.off()
##################################
# Effect size versus alpha basis #
##################################

pdf(width = 7, height = 12, file = "./visuals/effectsizeModeltype.pdf")
par(mfrow = c(5,2),
    mar = c(1,2,1,1),
    oma = c(8,5, 5, 1))
k <- 1

for(i in c(which(names(dat) == "truep_stanHMC_all.1"),
           which(names(dat) == "fpr_stan_all"),
           which(names(dat) == "truep_deseq_all.1"),
           which(names(dat) == "fpr_deseq_all"),   
           which(names(dat) ==  "truep_edgeR_all.1"),
           which(names(dat) == "fpr_edger_all"), 
           which(names(dat) ==  "truep_wilcoxBH_FDR_all.1"),
           which(names(dat) == "fpr_wilcox_all"), 
           which(names(dat) == "truep_ancom_all.1"),
           which(names(dat) == "fpr_ancom_all"))){
  if(k %in% c(1,3,5,7,9)){
  stripchart(dat[,i] ~ 
               dat$modelType + dat$effectsize,
             vertical = TRUE, 
             data = dat, 
             method = "jitter", 
             pch = 20, 
             col = #c(add.alpha(colz, alpha = 0.9),
                     add.alpha(colz, alpha = 0.5),
                    # add.alpha(colz, alpha = 0.2)),
             cex=1.5,
             xaxt="n",
             yaxt="n",
             ylab="",
             ylim=c(0,1),
             xlim = c(0,9.5),
             frame.plot = FALSE,
             xpd = NA
  )
  axis(side = 2, 
       at = c(0,0.5,1),
       labels = c(0,0.5,1),
       las = 2, 
       cex.axis = 1.8)
  
  boxplot(dat[,i] ~ 
            dat$modelType + dat$effectsize,
          las=2,
          yaxt="n",
          xaxt="n",
          outline=F,
          add=T,
          col =# c(add.alpha(colz, alpha = 0.9),
                  add.alpha(colz, alpha = 0.5),
                 # add.alpha(colz, alpha = 0.2)),
          xpd = NA,
          xlim = c(0,9.5),
          frame.plot = FALSE
  )
  
  abline(v = c(3.5,6.5),
         lty = 2)
  }else{
    stripchart(dat[,i] ~ 
                 dat$modelType,
               vertical = TRUE, 
               data = dat, 
               method = "jitter", 
               pch = 20, 
               col = #c(add.alpha(colz, alpha = 0.9),
                 add.alpha(colz, alpha = 0.5),
               # add.alpha(colz, alpha = 0.2)),
               cex=1.5,
               xaxt="n",
               yaxt="n",
               ylab="",
               ylim=c(0,1),
               xlim = c(0,3.5),
               frame.plot = FALSE,
               xpd = NA
    )
    axis(side = 2, 
         at = c(0,0.5,1),
         labels = c(0,0.5,1),
         las = 2, 
         cex.axis = 1.8)
    
    boxplot(dat[,i] ~ 
              dat$modelType,
            las=2,
            yaxt="n",
            xaxt="n",
            outline=F,
            add=T,
            col =# c(add.alpha(colz, alpha = 0.9),
              add.alpha(colz, alpha = 0.5),
            # add.alpha(colz, alpha = 0.2)),
            xpd = NA,
            xlim = c(0,3.5),
            frame.plot = FALSE
    )
    
  }
  if(k == 1){
    text(x = 4.5, 
       y = 1.3, 
       "True positive rate",
       xpd = NA, 
       cex = 1.8)
  
    text(x = -2.8, 
       y = -0.055,
       srt = 90,
       adj = 0.5,
       "Dirichlet-multinomial (HMC)",
       xpd = NA, 
       cex = 1.5, 
       pos = 4)
  }
  if(k == 2){
    text(x = 2, 
         y = 1.3, 
         "False positive rate",
         xpd = NA, 
         cex = 1.8)
  }
  if(k == 3){
    text(x = -2.8, 
         y = 0.27,
         srt = 90,
         adj = 0.5, 
         "DESeq2",
         xpd = NA, 
         cex = 1.5,
         pos = 4)
  }
  if(k == 5){
    text(x = -2.8, 
         y = 0.27,
         srt = 90,
         adj = 0.5, 
         "edgeR",
         xpd = NA, 
         cex = 1.5,
         pos = 4)
  }
  if(k == 7){
    text(x = -2.8, 
         y = 0.27,
         srt = 90,
         adj = 0.5, 
         "Wilcoxon",
         xpd = NA, 
         cex = 1.5,
         pos = 4)
  }
  if(k == 9){
    text(x = -2.8, 
         y = 0.27,
         srt = 90,
         adj = 0.5, 
         "ANCOM",
         xpd = NA, 
         cex = 1.5,
         pos = 4)
    
    text("Effect\nsize",
         x = -1,
         y = -0.2,
         xpd = NA,
         cex = 1.5)
    
    text("1.1",
         x = 2,
         y = -0.18,
         xpd = NA,
         cex = 1.5)
    text("1.5",
         x = 5,
         y = -0.18,
         xpd = NA,
         cex = 1.5)
    text("2",
         x = 8,
         y = -0.18,
         xpd = NA,
         cex = 1.5)
    
    #Legend
    rect(xleft = 4.5,
         xright = 5.5,
         ytop = -0.4,
         ybottom = -0.3,
         border = NA,
         col = add.alpha(colz[1], alpha = 0.9),
         xpd = NA)
    
    text("Pareto (0.7)",
         x = 5.7,
         y = -0.35,
         xpd = NA,
         pos = 4,
         cex = 1.5)
    
    rect(xleft = 9.5,
         xright = 10.5,
         ytop = -0.4,
         ybottom = -0.3,
         border = NA,
         col = add.alpha(colz[2], alpha = 0.9),
         xpd = NA)
    
    text("Pareto (4)",
         x = 10.7,
         y = -0.35,
         xpd = NA,
         pos = 4,
         cex = 1.5)
    
    rect(xleft = 14,
         xright = 15,
         ytop = -0.4,
         ybottom = -0.3,
         border = NA,
         col = add.alpha(colz[3], alpha = 0.9),
         xpd = NA)
    
    text("Equal",
         x = 15.2,
         y = -0.35,
         xpd = NA,
         pos = 4,
         cex = 1.5)
  }
  
  # if(k == 10){
  #   text("1.1",
  #        x = 2,
  #        y = -0.18,
  #        xpd = NA,
  #        cex = 1.5)
  #   text("1.5",
  #        x = 5,
  #        y = -0.18,
  #        xpd = NA,
  #        cex = 1.5)
  #   text("2",
  #        x = 8,
  #        y = -0.18,
  #        xpd = NA,
  #        cex = 1.5)
  # }
  # 
  k <- k + 1
}
dev.off()


pdf(width = 7, height = 12, file = "./visuals/effectsizeModeltypenew.pdf")
par(mfrow = c(5,2),
    mar = c(1,2,1,1),
    oma = c(8,5, 5, 1))
k <- 1

for(i in c(which(names(dat) == "truep_diff_all.1"),
           which(names(dat) == "fpr_diff_all"),
           which(names(dat) ==  "truep_aldexW_all.1"),
           which(names(dat) == "fpr_aldexW_all"), 
           which(names(dat) ==  "truep_aldexT_all.1"),
           which(names(dat) == "fpr_aldexT_all"), 
           which(names(dat) ==  "truep_stanVB_all.1"),
           which(names(dat) == "fpr_stanvb_all"), 
          
          which(names(dat) == "truep_jags_clr_all.1"),
          which(names(dat) == "fpr_jags_clr_all"))){
  if(k %in% c(1,3,5,7,9)){
    stripchart(dat[,i] ~ 
                 dat$modelType + dat$effectsize,
               vertical = TRUE, 
               data = dat, 
               method = "jitter", 
               pch = 20, 
               col = #c(add.alpha(colz, alpha = 0.9),
                 add.alpha(colz, alpha = 0.5),
               # add.alpha(colz, alpha = 0.2)),
               cex=1.5,
               xaxt="n",
               yaxt="n",
               ylab="",
               ylim=c(0,1),
               xlim = c(0,9.5),
               frame.plot = FALSE,
               xpd = NA
    )
    axis(side = 2, 
         at = c(0,0.5,1),
         labels = c(0,0.5,1),
         las = 2, 
         cex.axis = 1.8)
    
    boxplot(dat[,i] ~ 
              dat$modelType + dat$effectsize,
            las=2,
            yaxt="n",
            xaxt="n",
            outline=F,
            add=T,
            col =# c(add.alpha(colz, alpha = 0.9),
              add.alpha(colz, alpha = 0.5),
            # add.alpha(colz, alpha = 0.2)),
            xpd = NA,
            xlim = c(0,9.5),
            frame.plot = FALSE
    )
    
    abline(v = c(3.5,6.5),
           lty = 2)
  }else{
    stripchart(dat[,i] ~ 
                 dat$modelType,
               vertical = TRUE, 
               data = dat, 
               method = "jitter", 
               pch = 20, 
               col = #c(add.alpha(colz, alpha = 0.9),
                 add.alpha(colz, alpha = 0.5),
               # add.alpha(colz, alpha = 0.2)),
               cex=1.5,
               xaxt="n",
               yaxt="n",
               ylab="",
               ylim=c(0,1),
               xlim = c(0,3.5),
               frame.plot = FALSE,
               xpd = NA
    )
    axis(side = 2, 
         at = c(0,0.5,1),
         labels = c(0,0.5,1),
         las = 2, 
         cex.axis = 1.8)
    
    boxplot(dat[,i] ~ 
              dat$modelType,
            las=2,
            yaxt="n",
            xaxt="n",
            outline=F,
            add=T,
            col =# c(add.alpha(colz, alpha = 0.9),
              add.alpha(colz, alpha = 0.5),
            # add.alpha(colz, alpha = 0.2)),
            xpd = NA,
            xlim = c(0,3.5),
            frame.plot = FALSE
    )
    
  }
if(k == 1){
  text(x = 4.5, 
       y = 1.3, 
       "True positive rate",
       xpd = NA, 
       cex = 1.8)
  
  text(x = -2.8, 
       y = 0.27,
       srt = 90,
       adj = 0.5,
       "DMM (MCMC)",
       xpd = NA, 
       cex = 1.5, 
       pos = 4)
}
if(k == 2){
  text(x = 2, 
       y = 1.3, 
       "False positive rate",
       xpd = NA, 
       cex = 1.8)
}
if(k == 3){
  text(x = -2.8, 
       y = 0.17,
       srt = 90,
       adj = 0.5, 
       "ALDEx2 (Welch's t)",
       xpd = NA, 
       cex = 1.5,
       pos = 4)
}
if(k == 5){
  text(x = -2.8, 
       y = 0.17,
       srt = 90,
       adj = 0.5, 
       "ALDEx2 (Wilcoxon)",
       xpd = NA, 
       cex = 1.5,
       pos = 4)
}
if(k == 7){
  text(x = -2.8, 
       y = 0.27,
       srt = 90,
       adj = 0.5, 
       "DMM (VI)",
       xpd = NA, 
       cex = 1.5,
       pos = 4)
}
if(k == 9){
  text(x = -2.8, 
       y = 0.27,
       srt = 90,
       adj = 0.5, 
       "DMM (MCMC CLR)",
       xpd = NA, 
       cex = 1.5,
       pos = 4)
  
  text("Effect\nsize",
       x = -1,
       y = -0.2,
       xpd = NA,
       cex = 1.5)
  
  text("1.1",
       x = 2,
       y = -0.18,
       xpd = NA,
       cex = 1.5)
  text("1.5",
       x = 5,
       y = -0.18,
       xpd = NA,
       cex = 1.5)
  text("2",
       x = 8,
       y = -0.18,
       xpd = NA,
       cex = 1.5)
  
  #Legend
  rect(xleft = 4.5,
       xright = 5.5,
       ytop = -0.4,
       ybottom = -0.3,
       border = NA,
       col = add.alpha(colz[1], alpha = 0.9),
       xpd = NA)
  
  text("Pareto (0.7)",
       x = 5.7,
       y = -0.35,
       xpd = NA,
       pos = 4,
       cex = 1.5)
  
  rect(xleft = 9.5,
       xright = 10.5,
       ytop = -0.4,
       ybottom = -0.3,
       border = NA,
       col = add.alpha(colz[2], alpha = 0.9),
       xpd = NA)
  
  text("Pareto (4)",
       x = 10.7,
       y = -0.35,
       xpd = NA,
       pos = 4,
       cex = 1.5)
  
  rect(xleft = 14,
       xright = 15,
       ytop = -0.4,
       ybottom = -0.3,
       border = NA,
       col = add.alpha(colz[3], alpha = 0.9),
       xpd = NA)
  
  text("Equal",
       x = 15.2,
       y = -0.35,
       xpd = NA,
       pos = 4,
       cex = 1.5)
}
  k <- k + 1
}
dev.off()


# # ########
# # # Plot only rare or abund/etc.. ABANDONED
# # ########
# 
# dat2 <- dat[which(dat$modelType != "uniform"),]
# 
# pdf(width = 7, height = 12, file = "./visuals/effectsizeModeltypeRare.pdf")
# par(mfrow = c(5,2),
#     mar = c(1,2,1,1),
#     oma = c(8,5, 5, 1))
# k <- 1
# colz <- c("coral2",#RColorBrewer::brewer.pal(n = 5, name = "Accent")[1],
#           "gold3")
# 
#   for(i in c(which(names(dat2) == "truep_diff_rare.1"),
#            which(names(dat2) == "fpr_diff_rare"),
#            which(names(dat2) == "truep_deseq_rare.1"),
#            which(names(dat2) == "fpr_deseq_rare"),   
#            which(names(dat2) ==  "truep_edgeR_rare.1"),
#            which(names(dat2) == "fpr_edger_rare"), 
#            which(names(dat2) ==  "truep_wilcoxBH_FDR_rare.1"),
#            which(names(dat2) == "fpr_wilcox_rare"), 
#            which(names(dat2) == "truep_ancom_rare.1"),
#            which(names(dat2) == "fpr_ancom_rare"))){
#     if(k %in% c(1,3,5,7,9)){
#       stripchart(dat2[,i] ~
#                dat2$modelType + dat2$effectsize,
#              vertical = TRUE,
#              dat = dat2,
#              method = "jitter",
#              pch = 20,
#              col = #c(add.alpha(colz, alpha = 0.9),
#                add.alpha(colz, alpha = 0.5),
#              # add.alpha(colz, alpha = 0.2)),
#              cex=1.5,
#              xaxt="n",
#              yaxt="n",
#              ylab="",
#              ylim=c(0,1),
#              xlim = c(0,9.5),
#              frame.plot = FALSE,
#              xpd = NA
#   )
#   axis(side = 2,
#        at = c(0,0.5,1),
#        labels = c(0,0.5,1),
#        las = 2,
#        cex.axis = 1.8)
# 
#   boxplot(dat2[,i] ~
#             dat2$modelType + dat2$effectsize,
#           las=2,
#           yaxt="n",
#           xaxt="n",
#           outline=F,
#           add=T,
#           col =# c(add.alpha(colz, alpha = 0.9),
#             add.alpha(colz, alpha = 0.5),
#           # add.alpha(colz, alpha = 0.2)),
#           xpd = NA,
#           xlim = c(0,9.5),
#           frame.plot = FALSE
#   )
# 
#   abline(v = c(2.5,4.5),
#          lty = 2)
# }else{
#   stripchart(dat[,i] ~ 
#                dat$modelType,
#              vertical = TRUE, 
#              data = dat, 
#              method = "jitter", 
#              pch = 20, 
#              col = #c(add.alpha(colz, alpha = 0.9),
#                add.alpha(colz, alpha = 0.5),
#              # add.alpha(colz, alpha = 0.2)),
#              cex=1.5,
#              xaxt="n",
#              yaxt="n",
#              ylab="",
#              ylim=c(0,1),
#              xlim = c(0,3.5),
#              frame.plot = FALSE,
#              xpd = NA
#   )
#   axis(side = 2, 
#        at = c(0,0.5,1),
#        labels = c(0,0.5,1),
#        las = 2, 
#        cex.axis = 1.8)
#   
#   boxplot(dat[,i] ~ 
#             dat$modelType,
#           las=2,
#           yaxt="n",
#           xaxt="n",
#           outline=F,
#           add=T,
#           col =# c(add.alpha(colz, alpha = 0.9),
#             add.alpha(colz, alpha = 0.5),
#           # add.alpha(colz, alpha = 0.2)),
#           xpd = NA,
#           xlim = c(0,3.5),
#           frame.plot = FALSE
#   )
#   
# }
#   if(k == 1){
#     text(x = 4.5,
#          y = 1.3,
#          "% true positives \n recovered",
#          xpd = NA,
#          cex = 1.8)
# 
#     text(x = -2.8,
#          y = -0.055,
#          srt = 90,
#          adj = 0.5,
#          "Dirichlet-multinomial",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 2){
#     text(x = 2,
#          y = 1.3,
#          "False positive rate",
#          xpd = NA,
#          cex = 1.8)
#   }
#   if(k == 3){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "DESeq2",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 5){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "edgeR",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 7){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "Wilcoxon",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 9){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "ANCOM",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
# 
#     text("Effect\nsize",
#          x = -1,
#          y = -0.2,
#          xpd = NA,
#          cex = 1.5)
# 
#     text("1.1",
#          x = 1.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#     text("1.5",
#          x = 3.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#     text("2",
#          x = 5.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
# 
#     #Legend
#     rect(xleft = 4.5,
#          xright = 5.5,
#          ytop = -0.4,
#          ybottom = -0.3,
#          border = NA,
#          col = add.alpha(colz[1], alpha = 0.9),
#          xpd = NA)
# 
#     text("Pareto (0.7)",
#          x = 5.7,
#          y = -0.35,
#          xpd = NA,
#          pos = 4,
#          cex = 1.5)
# 
#     rect(xleft = 9.5,
#          xright = 10.5,
#          ytop = -0.4,
#          ybottom = -0.3,
#          border = NA,
#          col = add.alpha(colz[2], alpha = 0.9),
#          xpd = NA)
# 
#     text("Pareto (4)",
#          x = 10.7,
#          y = -0.35,
#          xpd = NA,
#          pos = 4,
#          cex = 1.5)
#   }
# 
#   # if(k == 10){
#   #   text("1.1",
#   #        x = 1.5,
#   #        y = -0.18,
#   #        xpd = NA,
#   #        cex = 1.5)
#   #   text("1.5",
#   #        x = 3.5,
#   #        y = -0.18,
#   #        xpd = NA,
#   #        cex = 1.5)
#   #   text("2",
#   #        x = 5.5,
#   #        y = -0.18,
#   #        xpd = NA,
#   #        cex = 1.5)
#   # }
# 
#   k <- k + 1
# }
# dev.off()
# 
# #Plot results for medium features
# 
# pdf(width = 7, height = 12, file = "./visuals/effectsizeModeltypeMed.pdf")
# par(mfrow = c(5,2),
#     mar = c(1,2,1,1),
#     oma = c(8,5, 5, 1))
# k <- 1
# colz <- c("coral2",#RColorBrewer::brewer.pal(n = 5, name = "Accent")[1], 
#           "gold3", 
#           "cyan4") #c("darkorange","darkorchid","cyan3")
# 
# for(i in c(which(names(dat) == "truep_diff_med.1"),
#            which(names(dat) == "fpr_diff_med"),
#            which(names(dat) == "truep_deseq_med.1"),
#            which(names(dat) == "fpr_deseq_med"),   
#            which(names(dat) ==  "truep_edgeR_med.1"),
#            which(names(dat) == "fpr_edger_med"), 
#            which(names(dat) ==  "truep_wilcoxBH_FDR_med.1"),
#            which(names(dat) == "fpr_wilcox_med"), 
#            which(names(dat) == "truep_ancom_med.1"),
#            which(names(dat) == "fpr_ancom_med"))){
#   stripchart(dat[,i] ~
#                dat$modelType + dat$effectsize,
#              vertical = TRUE,
#              data = dat,
#              method = "jitter",
#              pch = 20,
#              col = #c(add.alpha(colz, alpha = 0.9),
#                add.alpha(colz, alpha = 0.5),
#              # add.alpha(colz, alpha = 0.2)),
#              cex=1.5,
#              xaxt="n",
#              yaxt="n",
#              ylab="",
#              ylim=c(0,1),
#              xlim = c(0,9.5),
#              frame.plot = FALSE,
#              xpd = NA
#   )
#   axis(side = 2,
#        at = c(0,0.5,1),
#        labels = c(0,0.5,1),
#        las = 2,
#        cex.axis = 1.8)
#   
#   boxplot(dat[,i] ~
#             dat$modelType + dat$effectsize,
#           las=2,
#           yaxt="n",
#           xaxt="n",
#           outline=F,
#           add=T,
#           col =# c(add.alpha(colz, alpha = 0.9),
#             add.alpha(colz, alpha = 0.5),
#           # add.alpha(colz, alpha = 0.2)),
#           xpd = NA,
#           xlim = c(0,9.5),
#           frame.plot = FALSE
#   )
#   
#   abline(v = c(3.5,6.5),
#          lty = 2)
#   
#   if(k == 1){
#     text(x = 4.5,
#          y = 1.3,
#          "% true positives \n recovered",
#          xpd = NA,
#          cex = 1.8)
#     
#     text(x = -2.8,
#          y = -0.055,
#          srt = 90,
#          adj = 0.5,
#          "Dirichlet-multinomial",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 2){
#     text(x = 4.5,
#          y = 1.3,
#          "False positive rate",
#          xpd = NA,
#          cex = 1.8)
#   }
#   if(k == 3){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "DESeq2",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 5){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "edgeR",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 7){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "Wilcoxon",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#   }
#   if(k == 9){
#     text(x = -2.8,
#          y = 0.27,
#          srt = 90,
#          adj = 0.5,
#          "ANCOM",
#          xpd = NA,
#          cex = 1.5,
#          pos = 4)
#     
#     text("Effect\nsize",
#          x = -1,
#          y = -0.2,
#          xpd = NA,
#          cex = 1.5)
#     
#     text("1.1",
#          x = 1.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#     text("1.5",
#          x = 3.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#     text("2",
#          x = 5.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#     
#     #Legend
#     rect(xleft = 4.5,
#          xright = 5.5,
#          ytop = -0.4,
#          ybottom = -0.3,
#          border = NA,
#          col = add.alpha(colz[1], alpha = 0.9),
#          xpd = NA)
#     
#     text("Pareto (0.7)",
#          x = 5.7,
#          y = -0.35,
#          xpd = NA,
#          pos = 4,
#          cex = 1.5)
#     
#     rect(xleft = 9.5,
#          xright = 10.5,
#          ytop = -0.4,
#          ybottom = -0.3,
#          border = NA,
#          col = add.alpha(colz[2], alpha = 0.9),
#          xpd = NA)
#     
#     text("Pareto (4)",
#          x = 10.7,
#          y = -0.35,
#          xpd = NA,
#          pos = 4,
#          cex = 1.5)
#     
#     rect(xleft = 14,
#          xright = 15,
#          ytop = -0.4,
#          ybottom = -0.3,
#          border = NA,
#          col = add.alpha(colz[3], alpha = 0.9),
#          xpd = NA)
#     
#     text("Equal",
#          x = 15.2,
#          y = -0.35,
#          xpd = NA,
#          pos = 4,
#          cex = 1.5)
#   }
#   
#   if(k == 10){
#     text("1.1",
#          x = 1.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#     text("1.5",
#          x = 3.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#     text("2",
#          x = 5.5,
#          y = -0.18,
#          xpd = NA,
#          cex = 1.5)
#   }
#   
#   k <- k + 1
# }
# dev.off()
# 


