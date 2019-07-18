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
#write.csv(redo, file = "redo.csv", row.names = F)

#REMAKE the combo variable without simrep
dat <- dat[,-which(names(dat) == "combo")]

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

dat <- data.frame(dat, truepPropsAll, truepPropsabund, truepPropsmed, truepPropsrare, truepPropsmedhi)

#Combine model performance among reps by averaging, do this for all model types not just the DMM
out <- data.frame(matrix(ncol = 187, 
                         nrow = length(unique(dat$combo))))
k <- 1
for(i in unique(dat$combo)){
  d <- subset(dat, combo == i)
  out[k,] <- c(d[1, c(1:9)], apply(d[,10:187], 2, FUN = mean))
  k <- k + 1
}

names(out) <- names(dat)[c(1:9, 10:187)]

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

dim(out)
######################################################
# True pos comparison plots with marginal histograms #
######################################################

biplot.with.marginalhist <- function(x, 
                                     y, 
                                     xlab = "",
                                     ylab = "",
                                     lr,
                                     ur,
                                     ll,
                                     ul) {
  
  def.par <- par(no.readonly = TRUE)
  xhist <- hist(x, plot = FALSE, breaks = seq(0,1, by = 0.05))
  yhist <- hist(y, plot = FALSE, breaks = seq(0,1, by = 0.05))
  nf <-
    
    layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(3, 1), c(1, 3), TRUE)
  
  par(mar = c(7, 6, 1, 1),
      oma = c(1, 1, 1, 3))
  
  plot(
    x = x,
    y = y,
    xlab = "",
    ylab = "",
    #line = 5,
    axes = F,
    ylim = c(0, 1),
    xlim = c(0,1),
    pch = 21,
    col = add.alpha("cyan4", alpha = 0.9),
    bg = add.alpha("cyan4", alpha = 0.4),
    xpd = NA
  )
  
  rect(xleft = ul,
       ybottom = 0.65,
       xright = ur, #0.46,
       ytop = 0.95,
       col = add.alpha("cyan4",0.1),
       border = NA
  )
  
  compare <- y/x 
  compare[compare == "Inf"] <- 1.1 #recoding when dividing by zero to show numerator is bigger
  compare[compare == "NaN"] <- NA
  
  text(0.23,0.8,
       length(which(compare > 1)),
       col = add.alpha("cyan4",0.8),
       font = 2,
       cex = 3)
  
  rect(xleft = ll,#0.61,
       ybottom = 0.05,
       xright = lr, #0.95,
       ytop = 0.36,
       col = add.alpha("coral2",0.1),
       border = NA
  )
  
  text(0.8,0.2,
       length(which(compare < 1)),
       col = add.alpha("coral2",0.8),
       font = 2,
       cex = 3)
  
  mtext(side = 1, xlab, line = 3.5)
  mtext(side = 2, ylab, line = 3)
  
  axis(
    side = 2,
    at = formatC(
      seq(0,1, by = 0.2),
      #format = "e",
      digits = 2
    ),
    labels = formatC(
      seq(0,1, by = 0.2),
      # format = "e",
      digits = 2
    ),
    las = 2
  )
  
  axis(
    side = 1,
    line = 1,
    at = formatC(
      seq(0,1, by = 0.2),
      #format = "e",
      digits = 2
    ),
    labels = formatC(
      seq(0,1, by = 0.2),
      # format = "e",
      digits = 2
    )
  )
  abline(
    a = 0,
    b = 1,
    lty = 5,
    lwd = 2,
    col = "cyan4"
  )
  
  par(mar = c(0, 6, 1, 1))
  barplot(xhist$counts,
          axes = F,
          space = 0,
          xlim = c(0, 20),
          las = 2,
          col = "cyan4",
          border = "gray55")
  
  lines(x = c(0,20), 
        y = c(0,0),
        lty = 1,
        col = "black")
  axis(
    side = 2,
    cex = 2,
    labels = c("", "" , ""),
    at = c(0, max(xhist$counts)/2, max(xhist$counts)),
    xpd = NA,
    las = 3
  )
  text(
    x = c(-3, -3, -3),
    y = c(0, max(xhist$counts)/2, max(xhist$counts)),
    cex = 1,
    c(0, max(xhist$counts)/2, max(xhist$counts)),
    xpd = NA
    #srt = 270
  )
  
  par(mar = c(6.5, 0, 1, 1))
  barplot(
    yhist$counts,
    xaxt = "n",
    space = 0,
    horiz = TRUE,
    #las = 3,
    ylim = c(0, 20 ),
    col = "cyan4",
    border = "gray55"
  )
  lines(x = c(0,0), 
        y = c(0,20),
        lty = 1,
        col = "black")
  axis(
    side = 1,
    cex = 2,
    line = 0.5,
    labels = c("", "" , ""),
    at = c(0, max(yhist$counts) / 2, max(yhist$counts)),
    xpd = NA
  )
  text(
    y = c(-4.4, -4.4, -4.4),
    x = c(0, max(yhist$counts), max(yhist$counts)),
    cex = 1,
    c(0, max(yhist$counts), max(yhist$counts)),
    xpd = NA,
    srt = 0
  )
  
  par(def.par)
}


pdf(width = 5,
    height = 5,
    file = "./visuals/DMvsEdgeR_truepos.pdf")
biplot.with.marginalhist(x = dat$truep_edgeR_all.1,
                         y = dat$truep_stanHMC_all.1,
                         ylab = "Proportion true positives DMMM",
                         xlab = "Proportion true positives edgeR",
                         ul = 0,
                         ll = 0.58,
                         ur = 0.46,
                         lr = 1.1)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/DMvsWilcoxon_truepos.pdf")
biplot.with.marginalhist(x = dat$truep_wilcoxBH_FDR_all.1,
                         y = dat$truep_stanHMC_all.1,
                         ylab = "Proportion true positives DMM",
                         xlab = "Proportion true positives Wilcoxon",
                         ul = 0,
                         ll = 0.68,
                         ur = 0.46,
                         lr = 0.92)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/DMvsAncom_truepos.pdf")
biplot.with.marginalhist(dat$truep_ancom_all.1,
                         dat$truep_stanHMC_all.1,
                         ylab = "Proportion true positives DMM",
                         xlab = "Proportion true positives ANCOM",
                         ul = 0,
                         ll = 0.64,
                         ur = 0.46,
                         lr = 0.96)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/DMvsDeseq_truepos.pdf")
biplot.with.marginalhist(dat$truep_deseq_all.1,
                         dat$truep_stanHMC_all.1,
                         ylab = "Proportion true positives DMM",
                         xlab = "Proportion true positives DESeq2",
                         ul = 0,
                         ll = 0.61,
                         ur = 0.46,
                         lr = 0.98)
dev.off()

#comparison among ppd estimation methods

pdf(width = 5,
    height = 5,
    file = "./visuals/MCMCvsHMC_truepos.pdf")
biplot.with.marginalhist(dat$truep_diff_all.1,
                         dat$truep_stanHMC_all.1,
                         ylab = "Proportion true positives HMC",
                         xlab = "Proportion true positives MCMC",
                         ul = 0,
                         ll = 0.61,
                         ur = 0.46,
                         lr = 0.98)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/MCMCvsHMC_falsePos.pdf")
biplot.with.marginalhist(dat$fpr_diff_all,
                         dat$fpr_stan_all,
                         ylab = "Proportion false positives HMC",
                         xlab = "Proportion false positives MCMC",
                         ul = 0,
                         ll = 0.61,
                         ur = 0.46,
                         lr = 0.98)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/MCMCvsvb_falsePos.pdf")
biplot.with.marginalhist(dat$fpr_diff_all,
                         dat$fpr_stanvb_all,
                         ylab = "Proportion false positives VI",
                         xlab = "Proportion false positives MCMC",
                         ul = 0,
                         ll = 0.61,
                         ur = 0.46,
                         lr = 0.98)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/hmcvsvb_falsePos.pdf")
biplot.with.marginalhist(dat$fpr_stan_all,
                         dat$fpr_stanvb_all,
                         ylab = "Proportion false positives VI",
                         xlab = "Proportion false positives HMC",
                         ul = 0,
                         ll = 0.61,
                         ur = 0.46,
                         lr = 0.98)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/VIvsHMC_truepos.pdf")
biplot.with.marginalhist(dat$truep_stanVB_all.1,
                         dat$truep_stanHMC_all.1,
                         ylab = "Proportion true positives HMC",
                         xlab = "Proportion true positives VI",
                         ul = 0,
                         ll = 0.61,
                         ur = 0.46,
                         lr = 0.98)
dev.off()

pdf(width = 5,
    height = 5,
    file = "./visuals/VIvsMCMC_truepos.pdf")
biplot.with.marginalhist(
                         dat$truep_diff_all.1,
                         dat$truep_stanVB_all.1,
                         xlab = "Proportion true positives MCMC",
                         ylab = "Proportion true positives VI",
                         ul = 0,
                         ll = 0.61,
                         ur = 0.46,
                         lr = 0.98)
dev.off()

###########################################
#  boxplots showing effect of replication #
###########################################
# 
# out <- out[order(out$nreps),]
# pdf(width = 7, height = 10, file = "./visuals/replicationVsPerformance.pdf")
# 
# par(mfrow = c(5,2), 
#     mar = c(1,2,1,1),
#     oma = c(8,6, 5, 1))
# 
# stripchart(dat$truep_diff_all.1~ dat$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            ylim=c(0,1),
#            xlim = c(0,3.5),
#            frame.plot = FALSE,
#            xpd = NA
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2, 
#      cex.axis = 1.8)
# 
# boxplot(dat$truep_diff_all.1 ~ dat$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         xpd = NA,
#         xlim = c(0,3.5),
#         frame.plot = FALSE
# )
# 
# text(x = 2, 
#      y = 1.5, 
#      "% true positives \n recovered",
#      xpd = NA, 
#      cex = 1.8)
# 
# text(x = -1.2, 
#      y = -0.055,
#      srt = 90,
#      adj = 0.5,
#      "Dirichlet-multinomial",
#      xpd = NA, 
#      cex = 1.5, 
#      pos = 4)
# stripchart(out[,137] ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# 
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# 
# boxplot(out[,137] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# text(x = 2, 
#      y = 1.5, 
#      "False positive rate",
#      xpd = NA, 
#      cex = 2)
# 
# #deseq
# stripchart(out[,164]  ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# text(x = -1.2, 
#      y = 0.27,
#      srt = 90,
#      adj = 0.5, 
#      "DESeq2",
#      xpd = NA, 
#      cex = 1.5,
#      pos = 4)
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# 
# boxplot(out[,164]  ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# 
# stripchart(out[,147]  ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# 
# boxplot(out[,147] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# #edger
# stripchart(out[,165] ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# text(x = -1.2, 
#      y = 0.34,
#      srt = 90,
#      adj = 0.5, 
#      "edgeR",
#      xpd = NA, 
#      cex = 1.5,
#      pos = 4)
# boxplot(out[,165] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# stripchart(out[,152] ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# 
# boxplot(out[,152] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# #mannwhitney/wilcoxon
# stripchart(out[,163] ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# text(x = -1.2, 
#      y = 0.25,
#      srt = 90,
#      adj = 0.5,
#      "Wilcoxon",
#      xpd = NA, 
#      cex = 1.5,
#      pos = 4)
# boxplot(out[,163] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# stripchart(out[,142] ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# 
# boxplot(out[,142] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# #ancom
# stripchart(out[,166]  ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# text(x = -1.2, 
#      y = 0.27,
#      srt = 90,
#      adj = 0.5, 
#      "ANCOM",
#      xpd = NA, 
#      cex = 1.5,
#      pos = 4)
# boxplot(out[,166] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# text(cex=2, 
#      x = c(0.7,1.8),
#      y  = 0,
#      adj = c(1,2.4),
#      c(10,50), 
#      xpd=NA, 
#      srt=45)
# 
# stripchart(out[,157] ~ out$nreps,
#            vertical = TRUE, 
#            data = out, 
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha("gray", alpha = 0.7),
#            cex=2.5,
#            xaxt="n",
#            yaxt="n",
#            ylab="",
#            xlim = c(0,3.5),
#            ylim=c(0,1),
#            frame.plot = FALSE
# )
# axis(side = 2, 
#      at = c(0,0.5,1),
#      labels = c(0,0.5,1),
#      las = 2,
#      cex.axis = 1.8)
# 
# boxplot(out[,157] ~ out$nreps,
#         las=2,
#         yaxt="n",
#         xaxt="n",
#         outline=F,
#         add=T,
#         col = NA,
#         frame.plot = FALSE
# )
# 
# 
# text(cex=2, 
#      x = c(0.7,1.8),
#      y  = 0,
#      adj = c(1,2.4),
#      c(10,50), 
#      xpd=NA, 
#      srt=45)
# 
# mtext(side = 1, 
#       "Replicates",
#       line = 4.5,
#       adj = -0.8,
#       cex = 2)
# dev.off()

