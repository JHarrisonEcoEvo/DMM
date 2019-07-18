#!/usr/bin/Rscript

#J. Harrison
#December 2018 
#University of Wyoming

library(gtools)
library(rjags)
library(DESeq2)
library(edgeR)
library(EnvStats)
library(rstan)
library(exactRankTests)
library(nlme)
library(ggplot2)
library(ALDEx2)

# module load gcc/7.3.0
# module load r/3.5.0-py27
# module load r-rjags/4-6
# module load r-coda/0.19-1
# module load r-lattice/0.20-35
# module load r-rstan/2.17.2-py27
# module load r-rcpp
# module load r-rcppeigen

options(scipen=99)

#README

#You will have an additional file that has paramaters that will read one line at a time
#Each line will paramaterize a job, and all jobs will be run in concert.

#This program is made to be serially run via Slurm using the wrap_slurm_Moran_DMsim.pl script
#that script requires you to change a loop range so that the range matches 
#the lines of parameters in your input parameter file. 

#To those perusing this script, it is best to start at the bottom of the script 
#where the function "main" is fed options from STDIN. 
#Once you understand parameter input, then read function "main" with detours to 
#each of the many functions that "main" calls. 

#The parameters to be fed to this script are made via the 
#simParamMaker.R script, which is just an expand grid of 
#parameter vectors printed to an output file.

#####################
#Model specification#
#####################

community.model.level <- "model{
for(i in 1:N){
for(j in start[i]:end[i]){
datamatrix[j,] ~ dmulti(p[j,], nreads[j])
p[j,1:notus] ~ ddirch(pi[i,]*theta[i])
}

pi[i,1:notus] ~ ddirch(alpha)
theta[i] ~ dunif(0, 4000)
}

for(k in 1:notus){
alpha[k] <- 0.0000001
}
}"

#######################
#######################
#Function definitions##
#######################
#######################

##################################
#presented in alphabetical order #
##################################


aldexRun <- function(groupings, 
                     otus, 
                     fname = fname,
                     simcom_out = simcom_out){
  
  aldeMC <-  aldex(t(simcom_out$simulatedCommunity),
                   groupings, 
                   mc.samples=1000, 
                   test="t", 
                   effect=TRUE,     
                   include.sample.summary=FALSE, 
                   denom="all", 
                   verbose=FALSE)
  
  #t test
  positivesT <- which(aldeMC$we.eBH <= 0.05)
  if(length(positivesT) >= 1){
    negativesT <- otus[-which(otus %in% positivesT)]
  }else{
    negativesT <- otus
  }
  
  #Wilcoxon test
  positivesW <- which(aldeMC$wi.eBH <= 0.05)
  if(length(positivesW) >= 1){
    negativesW <- otus[-which(otus %in% positivesW)]
  }else{
    negativesW <- otus
  }
  
  abund_category <- simcom_out$abund_category
  different <- simcom_out$differentTaxa
  params <- simcom_out$Parameters_alpha1_2_pi1_2
  
  positives_diffT <- data.frame(positivesT, abund_category[positivesT])
  negatives_diffT <- data.frame(negativesT, abund_category[negativesT])
  
  # fp_effectsT <- params[1, positivesT[!(positivesT %in% different)]] - params[2, positivesT[!(positivesT %in% different)]]
  # 
  # write.csv(fp_effectsT, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/",     fname,       "_fp_effectsTaldex.csv"))
  
  positives_diffW <- data.frame(positivesW, abund_category[positivesW])
  negatives_diffW <- data.frame(negativesW, abund_category[negativesW])
  
  # fp_effectsW <- params[1, positivesW[!(positivesW %in% different)]] - params[2, positivesW[!(positivesW %in% different)]]
  # 
  # write.csv(fp_effectsW, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/", fname,  "_fp_effectsWaldex.csv"))
  
  return(
    list(
      performanceAllT = fpr(
        otus = otus,
        xpos = positivesT,
        xneg = negativesT,
        different = simcom_out$differentTaxa
      ),
      performance_abundT = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diffT$positives[which(positives_diffT$abund_category.positives. == "abund")],
        xneg = negatives_diffT$negatives[which(negatives_diffT$abund_category.negatives. == "abund")],
        different = simcom_out$exp_abund_differ
      ),
      performance_medT = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diffT$positives[which(positives_diffT$abund_category.positives. == "med")],
        xneg = negatives_diffT$negatives[which(negatives_diffT$abund_category.negatives. == "med")],
        different = simcom_out$exp_med_differ
      ),
      performance_rareT = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diffT$positives[which(positives_diffT$abund_category.positives. == "rare")],
        xneg = negatives_diffT$negatives[which(negatives_diffT$abund_category.negatives. == "rare")],
        different = simcom_out$exp_rare_differ
      ),
      performance_med_hiT = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diffT$positives[which(positives_diffT$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diffT$negatives[which(negatives_diffT$abund_category.negatives. %in% c("abund", "med"))],
        different = c(simcom_out$exp_med_differ, simcom_out$exp_abund_differ)
      ),
      #Wilcoxon
      performanceAllW = fpr(
        otus = otus,
        xpos = positivesW,
        xneg = negativesW,
        different = different
      ),
      performance_abundW = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diffW$positives[which(positives_diffW$abund_category.positives. == "abund")],
        xneg = negatives_diffW$negatives[which(negatives_diffW$abund_category.negatives. == "abund")],
        different = simcom_out$exp_abund_differ
      ),
      performance_medW = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diffW$positives[which(positives_diffW$abund_category.positives. == "med")],
        xneg = negatives_diffW$negatives[which(negatives_diffW$abund_category.negatives. == "med")],
        different = simcom_out$exp_med_differ
      ),
      performance_rareW = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diffW$positives[which(positives_diffW$abund_category.positives. == "rare")],
        xneg = negatives_diffW$negatives[which(negatives_diffW$abund_category.negatives. == "rare")],
        different = simcom_out$exp_rare_differ
      ),
      performance_med_hiW = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diffW$positives[which(positives_diffW$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diffW$negatives[which(negatives_diffW$abund_category.negatives. %in% c("abund", "med"))],
        different = c(simcom_out$exp_med_differ, simcom_out$exp_abund_differ)
      )
    )
  )
}

#################
#ANCOM function #
#################

# For debugging
# simcom_out <- simCom(
#   modeltype = modeltype,
#   notus = notus,
#   nreads = nreads,
#   nreps = nreps,
#   precision = precision,
#   seed = seed,
#   num_da = num_da,
#   effectsize = effectsize
# )

ancom <- function(simulatedData, 
                  groupMem, 
                  different, 
                  otus, 
                  fname = fname,
                  params = simcom_out$Parameters_alpha1_2_pi1_2,
                  abund_category,
                  rare_different = simcom_out$exp_rare_differ,
                  med_different = simcom_out$exp_med_differ,
                  abund_different = simcom_out$exp_abund_differ){
  #This function takes the simulated community, 
  #the vector defining group membership, 
  #the vector of taxa that differ between treatment groups,
  #and a vector from 1 to the total number of otus present within the data.
  
  simulatedData[simulatedData == 1] <- 0 
  #Ancom looks for zeros to determine what to analyze, so converting back to zero
  
  simulatedData2 <- data.frame(seq(1, length(groupMem), 1), simulatedData)
  names(simulatedData2)[1] <- "Sample.ID"
  
  metadat <- data.frame(simulatedData2$Sample.ID, groupMem)
  names(metadat)[1] <-  "Sample.ID"
  
  #Perform the ANCOM analysis
  ancom.out <- ANCOM.main(OTUdat = simulatedData2, 
                          Vardat = metadat,
                          main.var = "groupMem",
                          adjusted = F,
                          sig = 0.05, 
                          multcorr = 2, 
                          repeated = FALSE,
                          prev.cut = 0.99) 
  
  #Identify those OTUs that ANCOM found to differ between sampling groups,
  #and those OTUs that do not differ
  #Note the if statement accounts for the bug in ANCOM where it assigns significance to all or most
  #features. This bug is documented on the QIIME forums
  #See: https://forum.qiime2.org/t/ancom-low-w-taxa-identified-as-significant-
  #issues-workaround-ancom2-code-instructions/6040
  #It is unclear what is causing this bug as it is not stable. 
  #ANCOM does not fail for all data simulated with the same parameters.
  
  if(length(which(ancom.out$W.taxa$detected_0.9 == TRUE)) >= (dim(simulatedData)[2] * 0.9) ){
    positives <- which(ancom.out$W.taxa$W_stat > 0)
    negatives <- otus[-which(otus %in% positives)]
  }else{
    positives <- which(colnames(simulatedData) %in% ancom.out$W.taxa$otu.names[ancom.out$W.taxa$detected_0.9==TRUE])
    if(length(positives) >= 1){
      negatives <- otus[-which(otus %in% positives)]
    }else{
      negatives <- otus
    }
  }
  
  positives_diff <- data.frame(positives, abund_category[positives])
  negatives_diff <- data.frame(negatives, abund_category[negatives])
  
  # fp_effects <- params[1, positives[!(positives %in% different)]] - 
  #   params[2, positives[!(positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/",    fname, "_fp_effectsDeseq.csv"))
  return(
    list(
      performanceAll = fpr(
        otus = otus,
        xpos = positives,
        xneg = negatives,
        different = different
      ),
      performance_abund = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "abund")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "abund")],
        different = abund_different 
      ),
      performance_med = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "med")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "med")],
        different = med_different
      ),
      performance_rare = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "rare")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "rare")],
        different = rare_different
      ),
      performance_med_hi = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      )
    )
  )
}


#The first two functions are the ANCOM software, written by Mandal
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

#CLR
clr <- function(x){
  log(x) -  (1/length(x))*sum(log(x))
}

######################################################
#DESEQ, determine differential expression using DESeq2
######################################################
# simulatedData = simcom_out$simulatedCommunity
# groupMem = groupings
# parameters = parameters
# different = simcom_out$differentTaxa
# otus = otus
# abund_category = simcom_out$abund_category

deseq <- function(simulatedData, 
                  groupMem, 
                  parameters, 
                  params = simcom_out$Parameters_alpha1_2_pi1_2,
                  different, 
                  otus, 
                  simcom_out,
                  fname = fname,
                  abund_category,
                  rare_different = simcom_out$exp_rare_differ,
                  med_different = simcom_out$exp_med_differ,
                  abund_different = simcom_out$exp_abund_differ){
  #This function takes the simulated community, 
  #the vector defining group membership, 
  # the vector of parameters used to make the data,
  #the vector of taxa that differ between treatment groups (notzero), 
  #and a vector from 1 to the total number of otus present within the data.
  
  #build two vectors, one for the treatment condition, the other for "batch", 
  #which is a needed argument for deseq
  
  treatmentFactor <- data.frame(rep("batch1", dim(simulatedData)[1]), 
                                groupMem)
  
  names(treatmentFactor) <- c("batch", "condition")
  
  #Remove OTUs with few reads, as per recommendations
  #NOTE: we are not currently doing this to facilitate comparison among methods
  # newmat <- simulatedData[,-which(colSums(simulatedData) <= 100)]
  # simulatedData <- newmat
  
  
  #Perform DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = t(simulatedData),
                                colData = treatmentFactor,
                                design= ~condition)
  
  dds <- estimateSizeFactors(dds)
  
  #NOTE we use gene wise dispersion estimates because when replicates are used to estimate dispersion
  #the function sometimes fails, particularly when data are simulated from a uniform distribution. 
  #I think this is because there is not enough variation among replicates for those cases. 
  
  dds <- estimateDispersionsGeneEst(dds)
  
  dispersions(dds) <- mcols(dds)$dispGeneEst
  
  dds <- nbinomWaldTest(dds)
  
  dds_results <- results(dds)
  
  #IMPORTANT: As per the DESEq2 vignette guidelines for benchmarking, I set NA values to 1. 
  #The adjusted p values output by DESEq2 are set to NA if the software thinks an outlier is present that should
  #not be trusted, or if the mean number of sequences for that taxon is low enough that the results
  #cannot be trusted. This is of course subject to change by the user in arbitrary ways.
  
  ps <- ifelse(is.na(dds_results$padj), 1, dds_results$padj)
  
  #Multiple comparison correction defaults to a Benjamini-Hochberg FDR correction
  
  positives<-  which(ps <= 0.05)
  negatives <- which(ps > 0.05)
  
  positives_diff <- data.frame(positives, abund_category[positives])
  negatives_diff <- data.frame(negatives, abund_category[negatives])
  
  #extract effect sizes for false positives
  # fp_effects <- params[1, positives[!(positives %in% different)]] - params[2, positives[!(positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/",  fname,  "_fp_effectsDeseq.csv"))
  # # 
  return(
    list(
      performanceAll = fpr(
        otus = otus,
        xpos = positives,
        xneg = negatives,
        different = different
      ),
      performance_abund = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "abund")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "abund")],
        different = abund_different 
      ),
      performance_med = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "med")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "med")],
        different = med_different
      ),
      performance_rare = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "rare")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "rare")],
        different = rare_different
      ),
      performance_med_hi = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      ),
      positives_diff = positives_diff,
      negatives_diff = negatives_diff
    )
  )
}

#####################################################################################
#diffAbund: calculate differential relative abundance between control and treatment##
#####################################################################################

#For debugging
# SimulatedData = simcom_out$simulatedCommunity
# modelRun_out = modelRun_out$mcmcSamples
# groupings = groupings
# different = simcom_out$differentTaxa
# otus = otus
# uncommon = uncommon
# abund_category = simcom_out$abund_category
# exp_abund_differ = simcom_out$exp_abund_differ
# exp_med_differ = simcom_out$exp_med_differ
# exp_rare_differ = simcom_out$exp_rare_differ


diffAbund <- function(SimulatedData,
                      modelRun_out,
                      different,
                      groupings,
                      Parameters_alpha1_2_pi1_2,
                      otus,
                      fname = fname,
                      abund_category,
                      rare_different ,
                      med_different,
                      abund_different){
  
  #Otus is a sequence from 1 to the number of otus in simulated dataset
  
  calc_certain_diffsG <- function(mcmc_of_diffs){
    positives <- vector()
    negatives <- vector()
    
    for(i in 1:dim(mcmc_of_diffs)[1]){
      perc <- length(which(mcmc_of_diffs[i,] > 0)) / length(mcmc_of_diffs[i,])
      if(perc >= 0.95 | perc <= 0.05){
        positives <- c(positives, i)
      }else{
        negatives <- c(negatives, i)
      }
    }
    
    return(list(positives = positives, 
                negatives = negatives))
  }
  
  
  diffs <- modelRun_out$pi[1,,,1:2] - modelRun_out$pi[2,,,1:2]
  out <- calc_certain_diffsG(cbind(diffs[,,1], diffs[,,2]))
  
  #extract effect sizes for false positives
  # fp_effects <- Parameters_alpha1_2_pi1_2[1, out$positives[!(out$positives %in% different)]] - Parameters_alpha1_2_pi1_2[2, out$positives[!(out$positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/", fname, "_fp_effectsDMJags.csv", sep = ""))
  
  ##############################################
  #Determine performance by abundance category #
  ##############################################
  
  positives_diff <- data.frame(out$positives, abund_category[out$positives])
  negatives_diff <- data.frame(out$negatives, abund_category[out$negatives])
  
  #######################
  #apply CLR and do over
  #######################
  #concatenate so that mcmc iterations are different columns
  mcmc_clr_1 <- apply(cbind(modelRun_out$pi[1,,,1], modelRun_out$pi[1,,,2]),2,FUN = clr) #taxa as rows, mcmc samples as columns. iterating by column, so that we compute clr across taxa for each mcmc
  
  mcmc_clr_2 <- apply(cbind(modelRun_out$pi[2,,,1], modelRun_out$pi[2,,,2]),2,FUN = clr)
  
  diffs <- data.frame(matrix(nrow = dim(mcmc_clr_1)[1], ncol = dim(mcmc_clr_1)[2])) #taxa rows, mcmc samples columns, iterate by 
  
  for(i in 1:dim(mcmc_clr_1)[2]){
    diffs[,i] <- mcmc_clr_1[,i] - mcmc_clr_2[,i] #differences in mcmc samples for all taxa
  }
  outclr <- calc_certain_diffsG(diffs)
  
  #   #The above line indexes properly because otus always starts at 1
  
  #extract effect sizes for false positives
  # fp_effects <- Parameters_alpha1_2_pi1_2[1, outclr$positives[!(outclr$positives %in% different)]] - Parameters_alpha1_2_pi1_2[2, outclr$positives[!(outclr$positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/", fname, "_fp_effectsDMJagsCLR.csv", sep = ""))
  
  return(
    list(
      performanceAll = fpr(
        otus = otus,
        xpos = out$positives,
        xneg = out$negatives,
        different = different
      ),
      performance_abund = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "abund")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "abund")],
        different = abund_different 
      ),
      performance_med = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "med")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "med")],
        different = med_different
      ),
      performance_rare = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "rare")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "rare")],
        different = rare_different
      ),
      performance_med_hi = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      )
      ,positives_diff = positives_diff,
      performanceAllclr = fpr(
        otus = otus,
        xpos = outclr$positives,
        xneg = outclr$negatives,
        different = different
      )
    )
  )
}


diffStan <- function(SimulatedData,
                     different,
                     groupings,
                     Parameters_alpha1_2_pi1_2 = simcom_out$Parameters_alpha1_2_pi1_2,
                     otus,
                     fname = fname,
                     DM = DM,
                     abund_category = simcom_out$abund_category,
                     rare_different = simcom_out$exp_rare_differ,
                     med_different = simcom_out$exp_med_differ,
                     abund_different = simcom_out$exp_abund_differ){
  
  
  fitstan_NUTS <- rstan::sampling(DM, 
                                  data=list(datamatrix=SimulatedData, 
                                            nreps=nrow(SimulatedData), 
                                            notus=ncol(SimulatedData), 
                                            N=2, 
                                            start = c(min(which(groupings == unique(groupings)[1])),
                                                      min(which(groupings ==unique(groupings)[2]))),
                                            end = c(max(which(groupings == unique(groupings)[1])),
                                                    max(which(groupings == unique(groupings)[2])))), 
                                  chains=2, 
                                  control = list(max_treedepth = 10),
                                  warmup=1000, 
                                  iter=2000, 
                                  thin=2, 
                                  algorithm="NUTS", 
                                  cores=2, 
                                  pars<-c("pi","theta"),
                                  init =  initcalculatorStan(dat = SimulatedData, 
                                                             starts = c(min(which(groupings == unique(groupings)[1])),
                                                                        min(which(groupings == unique(groupings)[2]))), 
                                                             ends = c(max(which(groupings == unique(groupings)[1])),
                                                                      max(which(groupings == unique(groupings)[2])))),
                                  verbose=T) 
  
  est.pi<-rstan::extract(fitstan_NUTS,"pi")
  
  diffs_NUTS <- est.pi$pi[,1,] - est.pi$pi[,2,]
  
  calc_certain_diffs <- function(mcmc_of_diffs){
    positives <- vector()
    negatives <- vector()
    
    for(i in 1:dim(mcmc_of_diffs)[2]){ #note the different dimension compared to rjags
      perc <- length(which(mcmc_of_diffs[,i] > 0 ))/ length(mcmc_of_diffs[,i])      
      if(perc >= 0.95 | perc <= 0.05){
        positives <- c(positives, i)
      }else{
        negatives <- c(negatives, i)
      }
    }
    return(list(positives = positives, 
                negatives = negatives))
  }
  
  outNuts <- calc_certain_diffs(diffs_NUTS)
  
  #extract effect sizes for false positives
  # fp_effects <- Parameters_alpha1_2_pi1_2[1, outNuts$positives[!(outNuts$positives %in% different)]] - Parameters_alpha1_2_pi1_2[2, outNuts$positives[!(outNuts$positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/",fname,"_fp_effectsDMStan.csv", sep = ""))
  
  positives_diff <- data.frame(outNuts$positives, abund_category[outNuts$positives])
  negatives_diff <- data.frame(outNuts$negatives, abund_category[outNuts$negatives])
  
  #Calculate RMSE
  
  est.pi<-apply(rstan::extract(fitstan_NUTS,"pi")$pi,c(2,3),mean)
  true.pi<- Parameters_alpha1_2_pi1_2
  
  error.pi<-sqrt(mean((est.pi-true.pi)^2))
  
  est.pi <-  rstan::extract(fitstan_NUTS,"pi")
  
  ci1 <- apply(est.pi$pi[,1,], 2, FUN=HDIofMCMC)
  ci2 <- apply(est.pi$pi[,2,], 2, FUN=HDIofMCMC)
  
  successes <- 0
  fails <- 0
  
  #Color designation is used for plots, if those are desired.
  for(i in otus){
    if(Parameters_alpha1_2_pi1_2[1,i] > ci1[1,i] & Parameters_alpha1_2_pi1_2[1,i] < ci1[2,i]){
      successes <- successes + 1
    }else{
      fails <- fails + 1
    }
    if(Parameters_alpha1_2_pi1_2[2,i] > ci2[1,i] & Parameters_alpha1_2_pi1_2[2,i] < ci2[2,i]){  
      successes <- successes + 1
    }else{
      fails <- fails + 1
    }
  }
  
  perf <- successes/(successes + fails)
  
  write.csv(summary(fitstan_NUTS, 
                    pars = c("pi"), 
                    probs = c(0.025, 0.975))$summary,
            file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/", fname,  "_DiagnosticsDMStan.csv", sep = ""))
  
  return(
    list(
      thetahmc = mean(rstan::extract(fitstan_NUTS,"theta")$theta),
      performanceAll = fpr(
        otus = otus,
        xpos = outNuts$positives,
        xneg = outNuts$negatives,
        different = different
      ),
      performance_abund = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "abund")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "abund")],
        different = abund_different 
      ),
      performance_med = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "med")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "med")],
        different = med_different
      ),
      performance_rare = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "rare")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "rare")],
        different = rare_different
      ),
      performance_med_hi = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      ),
      positives_diff = positives_diff,
      negatives_diff = negatives_diff,
      rmse_hmc = error.pi,
      hdi_overlap = perf
    )
  )
}


diffStan_VB <- function(SimulatedData= simcom_out$simulatedCommunity,
                        different = simcom_out$differentTaxa,
                        groupings,
                        Parameters_alpha1_2_pi1_2 = simcom_out$Parameters_alpha1_2_pi1_2,
                        otus,
                        fname = fname,
                        DM = DM,
                        abund_category,
                        rare_different = simcom_out$exp_rare_differ,
                        med_different = simcom_out$exp_med_differ,
                        abund_different = simcom_out$exp_abund_differ){
  
  #This function originally defined in modelRun function
  inits <- initcalculatorStan(dat = SimulatedData,
                              starts = c(min(which(groupings == unique(groupings)[1])),
                                         min(which(groupings ==unique(groupings)[2]))),
                              ends = c(max(which(groupings == unique(groupings)[1])),
                                       max(which(groupings == unique(groupings)[2]))))
  fitstan_VB <-rstan::vb(DM, 
                         data=list("datamatrix"=SimulatedData, 
                                   "nreps"=nrow(SimulatedData), 
                                   "notus"=ncol(SimulatedData),
                                   "N"=2, 
                                   "start" = c(min(which(groupings == unique(groupings)[1])),
                                               min(which(groupings ==unique(groupings)[2]))),
                                   "end" = c(max(which(groupings == unique(groupings)[1])),
                                             max(which(groupings == unique(groupings)[2])))), 
                         algorithm="meanfield", 
                         output_samples=1000, 
                         check_data = T,
                         seed=123,
                         pars<-c("pi","theta")
                         #Note initialization seems to only make very modest speed improvements for vb
                         
                         ,init = list("p"= inits[[1]]$p,
                                      "pi"= inits[[1]]$pi
                                      # "alpha"= inits[[1]]$alpha,
                                      # "theta"= inits[[1]]$theta
                         )
  ) 
  
  est.pi<-extract(fitstan_VB,"pi")
  
  diffs_VB <- est.pi$pi[,1,] - est.pi$pi[,2,]
  
  calc_certain_diffs <- function(mcmc_of_diffs){
    positives <- vector()
    negatives <- vector()
    
    for(i in 1:dim(mcmc_of_diffs)[2]){
      perc <- length(which(mcmc_of_diffs[,i] > 0 ))/ length(mcmc_of_diffs[,i])      
      if(perc >= 0.95 | perc <= 0.05){
        positives <- c(positives, i)
      }else{
        negatives <- c(negatives, i)
      }
    }
    return(list(positives = positives, 
                negatives = negatives))
  }
  
  outVB <- calc_certain_diffs(diffs_VB)
  
  #extract effect sizes for false positives
  # fp_effects <- Parameters_alpha1_2_pi1_2[1, outVB$positives[!(outVB$positives %in% different)]] - 
  #   Parameters_alpha1_2_pi1_2[2, outVB$positives[!(outVB$positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/", fname,"_fp_effectsDMStanVB.csv", sep = ""))
  # 
  positives_diff <- data.frame(outVB$positives, abund_category[outVB$positives])
  negatives_diff <- data.frame(outVB$negatives, abund_category[outVB$negatives])
  
  #Calculate RMSE
  
  est.pi<-apply(est.pi$pi, c(2,3), mean)
  true.pi<-Parameters_alpha1_2_pi1_2
  
  error.pi<-sqrt(mean((est.pi-true.pi)^2))
  
  est.pi<-extract(fitstan_VB,"pi")
  
  #Truth in HDI
  ci1 <- apply(est.pi$pi[,1,], 2, FUN=HDIofMCMC)
  ci2 <- apply(est.pi$pi[,2,], 2, FUN=HDIofMCMC)
  
  successes <- 0
  fails <- 0
  
  #Color designation is used for plots, if those are desired.
  for(i in otus){
    if(Parameters_alpha1_2_pi1_2[1,i] > ci1[1,i] & Parameters_alpha1_2_pi1_2[1,i] < ci1[2,i]){
      successes <- successes + 1
    }else{
      fails <- fails + 1
    }
    if(Parameters_alpha1_2_pi1_2[2,i] > ci2[1,i] & Parameters_alpha1_2_pi1_2[2,i] < ci2[2,i]){  
      successes <- successes + 1
    }else{
      fails <- fails + 1
    }
  }
  
  perf <- successes/(successes + fails)
  
  return(
    list(
      thetavb =   mean(rstan::extract(fitstan_VB,"theta")$theta),
      performanceAll = fpr(
        otus = otus,
        xpos = c(na.omit(outVB$positives)),
        xneg = c(na.omit(outVB$negatives)),
        different = different
      ),
      performance_abund = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "abund")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "abund")],
        different = abund_different 
      ),
      performance_med = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "med")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "med")],
        different = med_different
      ),
      performance_rare = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "rare")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "rare")],
        different = rare_different
      ),
      performance_med_hi = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      ),
      positives_diff = positives_diff,
      negatives_diff = negatives_diff,
      rmse_vb = error.pi,
      hdi_overlap = perf
      
    )
  )
}
################
#edgeR function#
################

edgeRator <- function(simulatedData, 
                      groupings, 
                      params = simcom_out$Parameters_alpha1_2_pi1_2,
                      different, 
                      otus, 
                      fname = fname,
                      abund_category = simcom_out$abund_category,
                      rare_different = simcom_out$exp_rare_differ,
                      med_different = simcom_out$exp_med_differ,
                      abund_different = simcom_out$exp_abund_differ){
  
  #first convert to an object suitable for the edgeR functions (class DGEList)
  # 
  # newmat <- simulatedData[,-which(colSums(simulatedData) <= 100)]
  # simulatedData <- newmat
  
  #Convert to DGEList
  y <- DGEList(counts = t(simulatedData), 
               group = groupings)
  
  #Calculate normalization factors
  y <- calcNormFactors(y)
  
  #Estimate dispersion, but first make design matrix
  design <- model.matrix(~groupings)
  
  y <-  estimateCommonDisp(y, design)
  
  #Use exact test to determine differentially expressed elements
  et <- exactTest(y)
  
  fit <-  glmQLFit(y)
  tr <- glmQLFTest(fit)
  
  positives <- which(et$table[,3] <= 0.05)
  negatives <- which(et$table[,3] > 0.05)
  
  positives_diff <- data.frame(positives, abund_category[positives])
  negatives_diff <- data.frame(negatives, abund_category[negatives])
  
  #extract effect sizes for false positives
  # fp_effects <- params[1, positives[!(positives %in% different)]] - params[2, positives[!(positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/",fname, "_fp_effectsedgeR.csv", sep = ""))
  # 
  return(
    list(
      performanceAll = fpr(
        otus = otus,
        xpos = positives,
        xneg = negatives,
        different = different
      ),
      performance_abund = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "abund")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "abund")],
        different = abund_different 
      ),
      performance_med = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "med")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "med")],
        different = med_different
      ),
      performance_rare = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "rare")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "rare")],
        different = rare_different
      ),
      performance_med_hi = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      ),
      positives_diff = positives_diff,
      negatives_diff = negatives_diff
    )
  )
}

#Calculate truepositives, false positives, and true and false negatives.
fpr <- function(xpos, xneg, otus, different){
  #xpos is vector of indices of features that differed
  #xneg is vector of indices of features that did not differ
  return(list(
    truep = length(which(xpos %in% different)),
    falsep =
      length(xpos) - length(which(xpos %in% different)),
    truen = length(which(xneg %in% otus[-different])),
    falsen = length(which(different %in% xneg))
  ))
}


#Code from Kruschke's Doing Bayesian Data Analysis
HDIofMCMC = function(sampleVec, credMass=0.95) {
  
  # Computes highest density interval from a sample of representative values,
  
  # estimated as shortest credible interval.
  
  # Arguments:
  
  # sampleVec
  
  # is a vector of representative values from a probability distribution.
  
  # credMass
  
  # is a scalar between 0 and 1, indicating the mass within the credible
  
  # interval that is to be estimated.
  
  # Value:
  
  # HDIlim is a vector containing the limits of the HDI 
  sortedPts = sort(sampleVec)
  
  ciIdxInc = ceiling(credMass * length(sortedPts))
  
  nCIs = length(sortedPts) - ciIdxInc
  
  ciWidth = rep(0, nCIs) 
  
  for(i in 1:nCIs) {
    
    ciWidth[i]= sortedPts[i + ciIdxInc] - sortedPts[i]
    
  }
  
  HDImin = sortedPts[which.min(ciWidth)]
  
  HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
  
  HDIlim = c(HDImin, HDImax)
  
  return(HDIlim)
  
}


initcalculator <- function(dat, starts, ends){
  notusinit <- dim(dat)[2]
  
  #get initial values for p
  p <- array(dim = c(max(ends),
                     notusinit,
                     2)
  )
  for(j in 1:2){
    for(i in 1:max(ends)){
      p[i,,j] <- unlist(dat[i,1:notusinit] /sum(dat[i,1:notusinit]))
    }
  }
  
  #Generate initial values for alpha
  #Note that I add a small value to avoid infinite density errors.
  alphas <- 0.0000000001 +  apply(p[,1:notusinit,1], 2, mean)
  names(alphas) <- NULL
  alphas2 <- 0.0000000001 +  apply(p[,1:notusinit,2], 2, mean)
  names(alphas2) <- NULL
  
  #Generate initial values for pi
  pi <- array(dim = c(length(starts),
                      notusinit,
                      2))
  
  for(j in 1:2){
    for(i in 1:length(starts)){
      pi[i,,j] <- (apply(dat[starts[i]:ends[i],1:notusinit], 2, mean) + 1) /
        (sum(apply(dat[starts[i]:ends[i],1:notusinit], 2, mean) + 1))
      pi[i,,j] <- pi[i,,j] + 0.0000000001
    }
  }
  
  # #Generate initial values for theta
  # theta <- matrix(nrow = 2,
  #                 ncol = length(starts))
  # for(j in 1:2){
  #   for(i in 1:length(starts)){
  #     thetas <- alphas / pi[i,,1] #no reason to use both parts of array
  #     
  #     theta <- mean(thetas[which(thetas < 15000)])
  #     
  #     #Here we avoid initializing using high values of theta that seem 
  #     #unrealistic.
  #   }
  # }
  
  return(list(list("p" = p[,,1] ,
                   "pi" = pi[,,1]
                   # "alpha" =  c(alphas),
                   # "theta"  = theta + 0.0000000001
                   # ,"hypertheta" = theta + 0.0000000001
  ),
  list("p" = p[,,2] ,
       "pi" = pi[,,2]
       # "alpha" =  c(alphas2),
       #"theta"  = theta + 0.0000000001
       # ,"hypertheta" = theta + 0.0000000001
  )
  )
  )
} 


#Function to calculate initial values for MCMC chains
initcalculatorStan <- function(dat, starts, ends){
  notusinit <- dim(dat)[2]
  
  #get initial values for p
  p <- array(dim = c(max(ends),
                     notusinit,
                     2)
  )
  for(j in 1:2){
    for(i in 1:max(ends)){
      p[i,,j] <- unlist(dat[i,1:notusinit] /sum(dat[i,1:notusinit]))
    }
  }
  
  #Generate initial values for alpha
  #Note that I add a small value to avoid infinite density errors.
  alphas <- apply(p[,1:notusinit,1], 2, mean)
  names(alphas) <- NULL
  alphas2 <- apply(p[,1:notusinit,2], 2, mean)
  names(alphas2) <- NULL
  
  #Generate initial values for pi
  pi <- array(dim = c(length(starts),
                      notusinit,
                      2))
  
  for(j in 1:2){
    for(i in 1:length(starts)){
      pi[i,,j] <- (apply(dat[starts[i]:ends[i],1:notusinit], 2, mean)) /
        (sum(apply(dat[starts[i]:ends[i],1:notusinit], 2, mean)))
    }
  }
  
  # #Generate initial values for theta
  # theta <- matrix(nrow = 2,
  #                 ncol = length(starts))
  # for(j in 1:2){
  #   for(i in 1:length(starts)){
  #     thetas <- alphas / pi[i,,1] #There is no reason to use both parts of array
  #     
  #     theta <- mean(thetas[which(thetas < 15000)])
  #     
  #     #Here we avoid initializing using high values of theta that seem 
  #     #unrealistic.
  #   }
  # }
  
  return(list(list("p" = p[,,1] ,
                   "pi" = pi[,,1]
                   # "alpha" =  c(alphas),
                   #"theta"  = theta + 0.0000000001
                   # ,"hypertheta" = theta + 0.0000000001
  ),
  list("p" = p[,,2] ,
       "pi" = pi[,,2]
       # "alpha" =  c(alphas2),
       #"theta"  = theta + 0.0000000001
       # ,"hypertheta" = theta + 0.0000000001
  )
  )
  )
}
########################################################
########################################################
#                    Main function
#this function calls all the others (aside from runner)#
########################################################
########################################################

main <- function(parameters){
  
  print(parameters)
  
  #make matrix to hold results
  output <- matrix(0, ncol = 219, nrow = 2)
  
  ## write parameters to output file
  output[1, 1] <- "modelType"
  output[2, 1] <- parameters[1] #modelType element of {uniform, pareto_shape_0.7, pareto_shape_4}
  
  output[1, 2] <- "notus"
  output[2, 2] <- parameters[2] #notus
  
  output[1, 3]<- "nreads"
  output[2, 3] <- parameters[3] #nreads
  
  output[1, 4] <- "nreps"
  output[2, 4] <- parameters[4] #nreps
  
  output[1, 5] <- "precision_theta"
  output[2, 5] <- parameters[5] 
  
  output[1, 6] <- "num_taxa_different"
  output[2, 6] <- parameters[6]
  
  output[1, 7] <- "effectsize"
  output[2, 7] <- parameters[7] #effect size; difference in relative abundance for true positives between treatment group.
  
  #For convenient reference we assign parameter values to objects
  modeltype <- parameters[1]
  notus <- as.numeric(parameters[2])
  nreads <- as.numeric(parameters[3])
  nreps <- as.numeric(parameters[4])
  precision <- as.numeric(parameters[5])
  num_da <-  as.numeric(parameters[6])
  seed <- 666 #metal AF
  effectsize <- as.numeric(parameters[7])
  
  ###############
  #Simulate data#
  ###############
  
  #This function is exceptionally important, 
  #review it carefully as you peruse this script. 
  
  simcom_out <- simCom(
    modeltype = modeltype,
    notus = notus,
    nreads = nreads,
    nreps = nreps,
    precision = precision,
    seed = seed,
    num_da = num_da,
    effectsize = effectsize
  )
  
  otus <- seq(1, length(simcom_out$simulatedCommunity), 
              by = 1)
  
  #Make a vector for group membership
  groupings <- c(rep("group1", nreps / 2), 
                 rep("group2", nreps / 2))
  
  #Create a character vector to use in file names of output
  fname <- paste(
    modeltype,
    "number_otus_",
    notus,
    "reads_",
    nreads,
    "reps_",
    nreps,
    "theta_",
    precision,
    "effectsize_",
    effectsize,
    "_",
    parameters[8],
    ".csv",
    sep = ""
  )
  
  # write.csv(simcom_out$simulatedCommunity, paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/SimCom", fname, sep = ""))
  
  ##############
  # Model data #
  ##############
  ptm <- proc.time()
  
  #Yay for Gibbs
  modelRun_out <- modelRun(SimulatedData = simcom_out$simulatedCommunity, 
                           groupings = groupings, 
                           parameters = parameters,
                           simcom_out = simcom_out
  )

  dmmTime <- c(proc.time() - ptm)[3]
  
  output[1, 9] <-"number_iterAdapt"
  output[2, 9] <- modelRun_out$iter_needed
  
  #Compute model performance, not mcmc diagnostics
  #For mcmc-diagnostics see the modelRun function and the files it writes.
  
  #Note that a future implementation could calculate this without
  #suspect parameters included in calculations. 
  #However, average model rmse is generally excellent.
  #This is partially due to the very small proportion estimates for some taxa.
  
  perf <- modelPerf(estimated_pis = modelRun_out$mcmcSamples$pi, 
                    expected_pis = simcom_out$Parameters_alpha1_2_pi1_2[1:2,],
                    abund_cat = simcom_out$abund_category,
                    otus = otus)
  
  output[1, 10:11] <- c("RMSE","perc_overlap")
  output[2, 10:11] <- c(perf$rmseNZmean, perf$proportionOverlapping)
  
  ############################################
  # Perform differential expression analyses #
  ############################################
  
  #Determine differential expression using our method  #
  #This function is also very important.               #
  
  
  diffout <- diffAbund(SimulatedData = simcom_out$simulatedCommunity,
                       modelRun_out = modelRun_out$mcmcSamples,
                       groupings = groupings,
                       otus = otus,
                       Parameters_alpha1_2_pi1_2 = simcom_out$Parameters_alpha1_2_pi1_2,
                       different = simcom_out$differentTaxa,
                       fname = fname,
                       abund_category = simcom_out$abund_category,
                       rare_different = simcom_out$exp_rare_differ,
                       med_different = simcom_out$exp_med_differ,
                       abund_different = simcom_out$exp_abund_differ)
  
  #True/false positives, true/false negatives through difference method
  output[1, 12:15]  <- c("truep_diff_all",
                         "falsep_diff_all",
                         "truen_diff_all",
                         "falsen_diff_all")
  output[2, 12:15] <- c(unlist(diffout$performanceAll))
  
  #This is vestigial code from a time when I was interested in calculating
  #performance for different abundance classes. If this is desired, then
  #the function will need to be amended.
  
  # output[1, 16:19]  <- c("truep_diff_abund",
  #                        "falsep_diff_abund",
  #                        "truen_diff_abund",
  #                        "falsen_diff_abund")
  # output[2, 16:19] <- c(unlist(diffout$performance_abund))
  # 
  # output[1, 20:23]  <- c("truep_diff_med",
  #                        "falsep_diff_med",
  #                        "truen_diff_med",
  #                        "falsen_diff_med")
  # output[2, 20:23] <- c(unlist(diffout$performance_med))
  # 
  # output[1, 24:27]  <- c("truep_diff_rare",
  #                        "falsep_diff_rare",
  #                        "truen_diff_rare",
  #                        "falsen_diff_rare")
  # output[2, 24:27] <- c(unlist(diffout$performance_rare))
  
  #Determine differential expression using Frequentist methodology.
  #Warnings are suppressed because p values can't be computed for zero values
  #and consequently this function throws lots of warnings.
  ptm <- proc.time()
  print("test")
  wilcox_out <- suppressWarnings(wilcoxonTester(simulatedData = simcom_out$simulatedCommunity,
                                                group = groupings,
                                                otus = otus,
                                                fname = fname,
                                                Parameters_alpha1_2_pi1_2 = simcom_out$Parameters_alpha1_2_pi1_2,
                                                different = simcom_out$differentTaxa,
                                                abund_category = simcom_out$abund_category,
                                                rare_different = simcom_out$exp_rare_differ,
                                                med_different = simcom_out$exp_med_differ,
                                                abund_different = simcom_out$exp_abund_differ))
  wilcoxTime <- c(proc.time() - ptm)[3]
  output[1, 28:31] <- c("truep_wilcox",
                        "falsep_wilcox",
                        "truen_wilcox",
                        "falsen_wilcox")
  
  output[2, 28:31] <- unlist(wilcox_out$performanceAll)
  
  output[1, 32:35]  <- c("truep_wilcoxBH_FDR_all",
                         "falsep_wilcoxBH_FDR_all",
                         "truen_wilcoxBH_FDR_all",
                         "falsen_wilcoxBH_FDR_all")
  
  output[2, 32:35] <- unlist(wilcox_out$performanceBH)
  # 
  # output[1, 36:39] <- c("truep_wilcoxBH_FDR_abund",
  #                       "falsep_wilcoxBH_FDR_abund",
  #                       "truen_wilcoxBH_FDR_abund",
  #                       "falsen_wilcoxBH_FDR_abund")
  # 
  # output[2, 36:39] <- unlist(wilcox_out$performance_abundBH)
  # 
  # output[1, 40:43]  <- c("truep_wilcoxBH_FDR_med",
  #                        "falsep_wilcoxBH_FDR_med",
  #                        "truen_wilcoxBH_FDR_med",
  #                        "falsen_wilcoxBH_FDR_med")
  # 
  # output[2, 40:43] <- unlist(wilcox_out$performance_medBH)
  # 
  # output[1, 44:47]  <- c("truep_wilcoxBH_FDR_rare",
  #                        "falsep_wilcoxBH_FDR_rare",
  #                        "truen_wilcoxBH_FDR_rare",
  #                        "falsen_wilcoxBH_FDR_rare")
  # 
  # output[2, 44:47] <- unlist(wilcox_out$performance_rareBH)
  
  #Medhi for Wilcox BH is later in output file. 
  
  #Compute differential expression using DESeq2 
  ptm <- proc.time()
  deseqout <- deseq(simulatedData = simcom_out$simulatedCommunity, 
                    groupMem = groupings,
                    otus = otus,
                    fname = fname,
                    params = simcom_out$Parameters_alpha1_2_pi1_2, 
                    different = simcom_out$differentTaxa,
                    abund_category = simcom_out$abund_category,
                    rare_different = simcom_out$exp_rare_differ,
                    med_different = simcom_out$exp_med_differ,
                    abund_different = simcom_out$exp_abund_differ) 
  
  deseqTime <- c(proc.time() - ptm)[3]
  
  output[1, 48:51] <- c("truep_deseq_all", 
                        "falsep_deseq_all", 
                        "truen_deseq_all", 
                        "falsen_deseq_all")
  output[2, 48:51] <- unlist(deseqout$performanceAll)
  # 
  # output[1, 52:55] <- c("truep_deseq_abund", 
  #                       "falsep_deseq_abund", 
  #                       "truen_deseq_abund", 
  #                       "falsen_deseq_abund")
  # output[2, 52:55] <- unlist(deseqout$performance_abund)
  # 
  # output[1, 56:59] <- c("truep_deseq_med", 
  #                       "falsep_deseq_med", 
  #                       "truen_deseq_med", 
  #                       "falsen_deseq_med")
  # output[2, 56:59] <- unlist(deseqout$performance_med)
  # 
  # output[1, 60:63] <- c("truep_deseq_rare", 
  #                       "falsep_deseq_rare", 
  #                       "truen_deseq_rare", 
  #                       "falsen_deseq_rare")
  # output[2, 60:63] <- unlist(deseqout$performance_rare)
  
  #Compute differential expression using edgeR
  #Important note, the edgeR manual recommends removing low expressed genes, which we do
  ptm <- proc.time()
  edgerout <- edgeRator(simulatedData = simcom_out$simulatedCommunity,
                        groupings = groupings,
                        otus = otus,
                        fname = fname,
                        params = simcom_out$Parameters_alpha1_2_pi1_2, 
                        different = simcom_out$differentTaxa,
                        abund_category = simcom_out$abund_category,
                        rare_different = simcom_out$exp_rare_differ,
                        med_different = simcom_out$exp_med_differ,
                        abund_different = simcom_out$exp_abund_differ)
  edgerTime <- c(proc.time() - ptm)[3]
  
  output[1, 64:67]  <- c("truep_edgeR_all", 
                         "falsep_edgeR_all", 
                         "truen_edgeR_all", 
                         "falsen_edgeR_all")
  
  output[2, 64:67] <- unlist(edgerout$performanceAll)
  # 
  # output[1, 68:71]  <- c("truep_edgeR_abund", 
  #                        "falsep_edgeR_abund", 
  #                        "truen_edgeR_abund", 
  #                        "falsen_edgeR_abund")
  # 
  # output[2, 68:71] <- unlist(edgerout$performance_abund)
  # 
  # output[1, 72:75]  <- c("truep_edgeR_med", 
  #                        "falsep_edgeR_med", 
  #                        "truen_edgeR_med", 
  #                        "falsen_edgeR_med")
  # 
  # output[2, 72:75] <- unlist(edgerout$performance_med)
  # 
  # output[1, 76:79]  <- c("truep_edgeR_rare", 
  #                        "falsep_edgeR_rare", 
  #                        "truen_edgeR_rare", 
  #                        "falsen_edgeR_rare")
  # 
  # output[2, 76:79] <- unlist(edgerout$performance_rare)

  #Compute differential expression using ancom
  ptm <- proc.time()
  ancomout <- ancom(simulatedData = simcom_out$simulatedCommunity,
                    groupMem = groupings, 
                    otus = otus,
                    fname = fname,
                    params = simcom_out$Parameters_alpha1_2_pi1_2, 
                    different = simcom_out$differentTaxa,
                    abund_category = simcom_out$abund_category,
                    rare_different = simcom_out$exp_rare_differ,
                    med_different = simcom_out$exp_med_differ,
                    abund_different = simcom_out$exp_abund_differ) 
  ancomTime <- c(proc.time() - ptm)[3]
  
  output[1, 80:83] <- c("truep_ancom_all", 
                        "falsep_ancom_all", 
                        "truen_ancom_all", 
                        "falsen_ancom_all")
  
  output[2, 80:83] <- unlist(ancomout$performanceAll)
  # 
  # output[1, 84:87]  <- c("truep_ancom_abund", 
  #                        "falsep_ancom_abund", 
  #                        "truen_ancom_abund", 
  #                        "falsen_ancom_abund")
  # 
  # output[2, 84:87] <- unlist(ancomout$performance_abund)
  # 
  # output[1, 88:91]  <- c("truep_ancom_med", 
  #                        "falsep_ancom_med", 
  #                        "truen_ancom_med", 
  #                        "falsen_ancom_med")
  # 
  # output[2, 88:91] <- unlist(ancomout$performance_med)
  # 
  # output[1, 92:95]  <- c("truep_ancom_rare", 
  #                        "falsep_ancom_rare", 
  #                        "truen_ancom_rare", 
  #                        "falsen_ancom_rare")
  # 
  # output[2, 92:95] <- unlist(ancomout$performance_rare)
  
  print("About to start Stan.")
  # module load gcc/7.3.0
  # module load r/3.5.0-py27
  # module load r-rjags/4-6
  #Run Stan HMC model 
  DM <- rstan::stan_model("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/DM.stan", 
                          model_name="DM")
  
  ptm <- proc.time()
  stanHMCout <- diffStan(SimulatedData = simcom_out$simulatedCommunity,
                         different = simcom_out$differentTaxa,
                         groupings = groupings,
                         Parameters_alpha1_2_pi1_2 = simcom_out$Parameters_alpha1_2_pi1_2,
                         otus = otus,
                         DM = DM,
                         fname = fname,
                         abund_category = simcom_out$abund_category,
                         rare_different = simcom_out$exp_rare_differ,
                         med_different = simcom_out$exp_med_differ,
                         abund_different = simcom_out$exp_abund_differ)
  stanTime <- c(proc.time() - ptm)[3]
  print("finished hmc")
  
  
  output[1, 96:99] <- c("truep_stanHMC_all", 
                        "falsep_stanHMC_all", 
                        "truen_stanHMC_all", 
                        "falsen_stanHMC_all")
  
  output[2, 96:99] <- unlist(stanHMCout$performanceAll)
  # 
  # output[1, 100:103]  <- c("truep_stanHMC_abund", 
  #                          "falsep_stanHMC_abund", 
  #                          "truen_stanHMC_abund", 
  #                          "falsen_stanHMC_abund")
  # 
  # output[2, 100:103] <- unlist(stanHMCout$performance_abund)
  # 
  # output[1, 104:107]  <- c("truep_stanHMC_med", 
  #                          "falsep_stanHMC_med", 
  #                          "truen_stanHMC_med", 
  #                          "falsen_stanHMC_med")
  # 
  # output[2, 104:107] <- unlist(stanHMCout$performance_med)
  # 
  # output[1, 108:111]  <- c("truep_stanHMC_rare", 
  #                          "falsep_stanHMC_rare", 
  #                          "truen_stanHMC_rare", 
  #                          "falsen_stanHMC_rare")
  # 
  # output[2, 108:111] <- unlist(stanHMCout$performance_rare)
  # 
  output[1, 112] <- "rmse_hmc"
  output[2, 112] <- stanHMCout$rmse_hmc
  
  ptm <- proc.time()
  stanVBout <- diffStan_VB(SimulatedData = simcom_out$simulatedCommunity,
                           different = simcom_out$differentTaxa,
                           groupings,
                           Parameters_alpha1_2_pi1_2 = simcom_out$Parameters_alpha1_2_pi1_2,
                           otus,
                           DM = DM,
                           fname = fname,
                           abund_category = simcom_out$abund_category,
                           rare_different = simcom_out$exp_rare_differ,
                           med_different = simcom_out$exp_med_differ,
                           abund_different = simcom_out$exp_abund_differ)
  
  stanVBTime <- c(proc.time() - ptm)[3]
  print("finished vb")
  output[1, 113:116] <- c("truep_stanVB_all", 
                          "falsep_stanVB_all", 
                          "truen_stanVB_all", 
                          "falsen_stanVB_all")
  
  output[2, 113:116] <- unlist(stanVBout$performanceAll)
  # 
  # output[1, 117:120]  <- c("truep_stanVB_abund", 
  #                          "falsep_stanVB_abund", 
  #                          "truen_stanVB_abund", 
  #                          "falsen_stanVB_abund")
  # 
  # output[2, 117:120] <- unlist(stanVBout$performance_abund)
  # 
  # output[1, 121:124]  <- c("truep_stanVB_med", 
  #                          "falsep_stanVB_med", 
  #                          "truen_stanVB_med", 
  #                          "falsen_stanVB_med")
  # 
  # output[2, 121:124] <- unlist(stanVBout$performance_med)
  # 
  # output[1, 125:128]  <- c("truep_stanVB_rare", 
  #                          "falsep_stanVB_rare", 
  #                          "truen_stanVB_rare", 
  #                          "falsen_stanVB_rare")
  # 
  # output[2, 125:128] <- unlist(stanVBout$performance_rare)  
  
  output[1, 129] <- "rmse_vb"
  output[2, 129] <- stanVBout$rmse_vb
  #} #Closes else statement started way back in the day: Line 460
  
  #Note these samples of theta do not approximate the theta used 
  #to simulate data. These represent noise within the data, that is, among replicate variation.
  
  output[1, 130] <- "thetaJAGS"
  output[2, 130] <- mean(apply(modelRun_out$mcmcSamples$theta[1,,1:2], 1, mean))
  
  output[1, 131] <- "actualOTUs"
  output[2, 131] <- dim(simcom_out$simulatedCommunity)[[2]]
  
  output[1, 132] <- "numDiffering"
  output[2, 132] <- length(simcom_out$differentTaxa)
  
  output[1, 133] <- "numDiffering_abund"
  output[2, 133] <- length(simcom_out$exp_abund_differ)
  
  output[1, 134] <- "numDiffering_med"
  output[2, 134] <- length(simcom_out$exp_med_differ)
  
  output[1, 135] <- "numDiffering_rare"
  output[2, 135] <- length(simcom_out$exp_rare_differ)
  
  output[1, 136] <- "simRep"
  output[2, 136] <- parameters[8] #no longer used
  
  output[1, 137] <- "rmseAbund"
  output[2, 137] <- perf$rmseAbund
  
  output[1, 138] <- "rmseMed"
  output[2, 138] <- perf$rmseMed
  
  output[1, 139] <- "rmseRare"
  output[2, 139] <- perf$rmseRare
  # 
  # #Late additon, combo of both medium and abundant results
  # output[1, 140:143]  <- c("truep_diff_medhi", 
  #                          "falsep_diff_medhi", 
  #                          "truen_diff_medhi", 
  #                          "falsen_diff_medhi")
  # output[2, 140:143] <- c(unlist(diffout$performance_med_hi))
  # 
  # output[1, 144:147]  <- c("truep_ancom_medhi", 
  #                          "falsep_ancom_medhi", 
  #                          "truen_ancom_medhi", 
  #                          "falsen_ancom_medhi")
  # output[2, 144:147] <- c(unlist(ancomout$performance_med_hi))
  # 
  # output[1, 148:151]  <- c("truep_wilcox_medhi", 
  #                          "falsep_wilcox_medhi", 
  #                          "truen_wilcox_medhi", 
  #                          "falsen_wilcox_medhi")
  # output[2, 148:151] <- c(unlist(wilcox_out$performance_med_hi))
  # 
  # output[1, 152:155]  <- c("truep_deseq_medhi", 
  #                          "falsep_deseq_medhi", 
  #                          "truen_deseq_medhi", 
  #                          "falsen_deseq_medhi")
  # output[2,  152:155] <- c(unlist(deseqout$performance_med_hi))
  # 
  # output[1, 156:159]  <- c("truep_edger_medhi", 
  #                          "falsep_edger_medhi", 
  #                          "truen_edger_medhi", 
  #                          "falsen_edger_medhi")
  # output[2,  156:159] <- c(unlist(edgerout$performance_med_hi))
  # 
  # output[1, 160:163]  <- c("truep_wilcoxBH_FDR_medhi", 
  #                          "falsep_wilcoxBH_FDR_medhi", 
  #                          "truen_wilcoxBH_FDR_medhi", 
  #                          "falsen_wilcoxBH_FDR_medhi")
  
  # output[2, 160:163] <- unlist(wilcox_out$performance_med_hiBH)
  
  output[1, 164]  <- "dmm_time"
  output[2,  164] <- dmmTime
  
  output[1, 165]  <- "wilcox_time"
  output[2,  165] <- wilcoxTime
  
  output[1, 166]  <- "deseq_time"
  output[2,  166] <- deseqTime
  
  output[1, 167]  <- "edger_time"
  output[2,  167] <- edgerTime
  
  output[1, 168]  <- "ancom_time"
  output[2,  168] <- ancomTime
  
  output[1, 169]  <- "stanHMCtime"
  output[2,  169] <- stanTime
  
  output[1, 170]  <- "stanVBtime"
  output[2,  170] <- stanVBTime
  
  output[1, 171:174]  <- c("truep_jags_clr_all", 
                           "falsep_jags_clr_all", 
                           "truen_jags_clr_all", 
                           "falsen_jags_clr_all")
  output[2, 171:174]  <- c(unlist(diffout$performanceAllclr))
  
  output[1, 175] <- "thetaHMC"
  output[2, 175] <- stanHMCout$thetahmc
  
  output[1, 176] <- "thetaVB_1"
  output[2, 176] <- stanVBout$thetavb
  
  ptm <- proc.time()
  print("aldex is about to start")
  aldexOut <- aldexRun(groupings= groupings, 
                       otus = otus, 
                       fname = fname,
                       simcom_out = simcom_out)
  
  aldexTime <- c(proc.time() - ptm)[3]
  
  output[1, 177:180] <- c("truep_aldexT_all", 
                          "falsep_aldexT_all", 
                          "truen_aldexT_all", 
                          "falsen_aldexT_all")
  
  output[2, 177:180] <- unlist(aldexOut$performanceAllT)
  # 
  # output[1, 181:184] <- c("truep_aldexT_abund", 
  #                         "falsep_aldexT_abund", 
  #                         "truen_aldexT_abund", 
  #                         "falsen_aldexT_abund")
  # 
  # output[2, 181:184] <- unlist(aldexOut$performance_abundT)
  # 
  # output[1, 185:188] <- c("truep_aldexT_med", 
  #                         "falsep_aldexT_med", 
  #                         "truen_aldexT_med", 
  #                         "falsen_aldexT_med")
  # 
  # output[2, 185:188] <- unlist(aldexOut$performance_medT)
  # 
  # output[1, 189:192] <- c("truep_aldexT_rare", 
  #                         "falsep_aldexT_rare", 
  #                         "truen_aldexT_rare", 
  #                         "falsen_aldexT_rare")
  # 
  # output[2, 189:192] <- unlist(aldexOut$performance_rareT)
  # 
  # output[1, 193:196] <- c("truep_aldexT_med_hi", 
  #                         "falsep_aldexT_med_hi", 
  #                         "truen_aldexT_med_hi", 
  #                         "falsen_aldexT_med_hi")
  # 
  # output[2, 193:196] <- unlist(aldexOut$performance_med_hiT)
  # 
  # output[1, 197:200] <- c("truep_aldexW_all", 
  #                         "falsep_aldexW_all", 
  #                         "truen_aldexW_all", 
  #                         "falsen_aldexW_all")
  # 
  # output[2,  197:200] <- unlist(aldexOut$performanceAllW)
  # 
  # output[1, 201:204] <- c("truep_aldexW_abund", 
  #                         "falsep_aldexW_abund", 
  #                         "truen_aldexW_abund", 
  #                         "falsen_aldexW_abund")
  # 
  # output[2, 201:204] <- unlist(aldexOut$performance_abundW)
  # 
  # output[1, 205:208] <- c("truep_aldexW_med", 
  #                         "falsep_aldexW_med", 
  #                         "truen_aldexW_med", 
  #                         "falsen_aldexW_med")
  # 
  # output[2, 205:208] <- unlist(aldexOut$performance_medW)
  # 
  # output[1, 209:212] <- c("truep_aldexW_rare", 
  #                         "falsep_aldexW_rare", 
  #                         "truen_aldexW_rare", 
  #                         "falsen_aldexW_rare")
  # 
  # output[2, 209:212] <- unlist(aldexOut$performance_rareW)
  # 
  # output[1, 213:216] <- c("truep_aldexW_med_hi", 
  #                         "falsep_aldexW_med_hi", 
  #                         "truen_aldexW_med_hi", 
  #                         "falsen_aldexW_med_hi")
  # 
  # output[2, 213:216] <- c(unlist(aldexOut$performance_med_hiW))
  
  output[1, 217] <- "aldexTime"
  output[2, 217] <- aldexTime
  
  
  output[1, 218] <- "truth_in_hdi_hmc"
  output[2, 218] <- stanHMCout$hdi_overlap
  
  output[1, 219] <- "truth_in_hdi_vb"
  output[2, 219] <- stanVBout$hdi_overlap
  
  print("yay done")
  ##write output
  write.table(output,file = paste( "/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/v4", fname,  sep = ""), col.names=F, row.names = F )
}

#############################
#calculate model performance#
#############################

modelPerf <- function(estimated_pis,
                      expected_pis,
                      abund_cat, 
                      otus = otus){
  #estimated_pis is the estimated pi from the model
  #expected_pis is the expected pi that was used to generate the simulated data
  #features_expected_to_differ are the features that were expected to differ between treatment groups
  #abund_cat describes which abundance class each feature belongs to. 
  
  #IMPORTANT NOTE: I do not remove rare (<100 observation) features here. 
  #It may be worth considering doing that, because RMSE may be higher for more abundant features
  
  #calculate model rmse
  rmse <- function(obs, exp){
    sqrt(sum((obs - exp)^2) / length(obs))
  }
  
  #Calculate rmse for both treatment groups.
  #Average across treatment groups, and get confidence intervals
  #not computing confidence intervals for now
  
  rmseR <- function(x){
    rmseAbund1 <- rmse(apply(estimated_pis[1, x,,1:2], 
                             1, 
                             mean),
                       expected_pis[1,x])
    rmseAbund2 <- rmse(apply(estimated_pis[2, x,,1:2], 
                             1, 
                             mean),
                       expected_pis[2,x])
    
    rmseNZmean <- mean(c(rmseAbund1, rmseAbund2))
    return(rmseNZmean)
  }
  
  ###############################################################################
  #Calculate model performance based on confidence interval overlap of 1:1 line##
  ###############################################################################
  
  #If the CI overlaps the expected value then we deem this evidence of successful
  #modeling. 
  
  ci1 <- apply(estimated_pis[1,,,1:2], 1, HDIofMCMC, credMass = 0.9)
  ci2 <- apply(estimated_pis[2,,,1:2], 1, HDIofMCMC, credMass = 0.9)
  
  #Build some matrices to hold output of following loop
  unpredicted <- matrix("black",
                        nrow=dim(expected_pis)[1], 
                        ncol=dim(expected_pis)[2])
  
  n_otus <- dim(expected_pis)[2]
  successes <- 0
  fails <- 0
  
  for(i in 1:n_otus){
    if(expected_pis[1,i] > ci1[1,i] & expected_pis[1,i] < ci1[2,i]){
      successes <- successes + 1
    }else{
      fails <- fails + 1
      unpredicted[1,i] <- "red"
    }
    if(expected_pis[2,i] > ci2[1,i] & expected_pis[2,i] < ci2[2,i]){  
      successes <- successes + 1
    }else{
      fails <- fails + 1
      unpredicted[2,i] <- "red"
    }
  }
  
  perf <- successes/(successes + fails)
  
  return(list(rmseNZmean = rmseR(otus),
              proportionOverlapping = perf, 
              unpredicted = unpredicted,
              rmseAbund = rmseR(grep("abund",abund_cat)),
              rmseMed = rmseR(grep("med",abund_cat)),
              rmseRare = rmseR(grep("rare",abund_cat))
  ))
}        

#Compute the Gelman-Rubin statistic
#Get Gelman-Rubin values for pis, keep track of the suspect taxa.
#Calculate MCMC diagnostics
mcmcdiag <- function(x, nparams) {
  #x is an mcmc object
  #nparams is number of params in the object
  Gr <- vector(length = nparams)
  GK <- vector(length = nparams)
  k <- 1
  a <- character(0)
  
  while (k <= nparams) {
    m <- x[1:length(x)][, k]
    
    gr <- gelman.diag(m) 
    print(paste("Feature", k, sep = " "))
    print("Gelman-Rubin")
    print(gr)
    if (gr[[1]][1] <= 2) {
      Gr[k] <- "passed"
    } else{
      Gr[k] <- "failed"
    }
    
    gk <- geweke.diag(m,
                      frac1 = 0.1,
                      frac2 = 0.5)
    suspectGK <-  names(which(2 * pnorm(-abs(gk[[1]]$z)) < 0.08))
    if (identical(a, suspectGK)) {
      GK[k] <- "passed"
    } else if (suspectGK == "var1") {
      GK[k] <- "failed"
    }
    
    k <- k + 1
  }
  return(list(Gr,
              GK))
}

#################################
#run Dirichlet-multinomial model#
#################################

#Functions are defined in place inside this function, for readability.
modelRun <- function(SimulatedData = simcom_out$simulatedCommunity, 
                     groupings, 
                     parameters, 
                     fname = fname,
                     simcom_out 
                    ){
  
  #This function takes the outputs from simCom and the parameters for the simulated data
  #groupings is the membership of each sample to a sampling group
  #parameters are the parameters used to simulate data
  #simcom_out is the output of the simulation function
  
  #compile model
  sim.mod.jags <- jags.model(
    textConnection(community.model.level),
    data = list(
      datamatrix = SimulatedData,
      notus = dim(SimulatedData)[2],
      nreads = rowSums(SimulatedData),
      N = length(unique(groupings)),
      start = c(min(which(
        groupings == unique(groupings)[1])),
        min(which(
          groupings == unique(groupings)[2]))),
      end = c(max(which(
        groupings == unique(groupings)[1])),
        max(which(
          groupings == unique(groupings)[2])))
    ),
    inits = initcalculator(dat = SimulatedData,
                           starts = c(min(which(groupings == unique(groupings)[1])),
                                      min(which(groupings == unique(groupings)[2]))), 
                           
                           ends = c(max(which(groupings == unique(groupings)[1])),
                                    max(which(groupings == unique(groupings)[2])))),
    n.chains = 2,
    n.adapt = 0
  )
  
  #Adapt model for as long as needed and save the adaptation period
  #This was required to avoid some models needing vastly longer adaptation periods than others.
  #Recall that adapt is a way to determine how to stride through the parameter space and is 
  #different than burnin
  # 
  iter_needed <- 0
  y = FALSE
  while(y == FALSE){
    y <-  adapt(sim.mod.jags,
                n.iter = 1000,
                end.adaptation = FALSE)
    iter_needed <- 1000 + iter_needed
    if(iter_needed > 20000){break}
  }
  
  #burn in model
  update(sim.mod.jags,
         n.iter = 300000)
  
  #here we extract the MCMC samples 
  sim.mod.sam <- jags.samples(model = sim.mod.jags, 
                              variable.names = c(
                                "pi", 
                                "theta",
                                "p"
                              ),
                              n.iter = 2000,
                              thin = 4)
  
  zed <- as.mcmc.list(sim.mod.sam$pi)
  
  #What have we wrought ??
  diagout <- mcmcdiag(zed, dim(SimulatedData)[2])
  suspect <- rep("failed", dim(SimulatedData)[2])
  suspect[which(diagout[[1]] == "passed")] <- "passed"
  
  suspectGK <- rep("failed", dim(SimulatedData)[2])
  suspectGK[which(diagout[[2]] == "passed")] <- "passed"
  
  mcmc_diag <- rbind(
    suspect,
    suspectGK
  )
  colnames(mcmc_diag) <- paste("feature", seq(1, dim(SimulatedData)[2], by = 1))
  
  # write MCMC diagnostics to disk
  ##write.csv(
  #   mcmc_diag,
  #   file = paste(
  #     "/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/MCMCdiag",
  #     fname,
  #     sep = ""
  #   ),
  #   row.names = F
  # )
  
  return(list(mcmcSamples = sim.mod.sam,
              iter_needed = iter_needed,
              mcmc_diag = mcmc_diag
  ))
}

#########
#Plotter#
#########

plotter <- function(estimatedVals, expectedVals, plotcolors){
  
  ws.q<-apply(estimatedVals[1,,,1:2], 1, quantile, probs = c(0.025, 0.975))
  
  plot(expectedVals[1,], apply(estimatedVals[1,,,1:2], 1, mean),
       main="",
       xlab="Expected pi",
       ylab="Estimated pi",
       ylim=c(min(ws.q[1,]), max(ws.q[2,])),
       col=plotcolors[1,])
  
  segments(expectedVals[1,], 
           ws.q[1,], 
           expectedVals[1,],
           ws.q[2,], 
           col=plotcolors[1,])
  
  #add other treatment group
  points(expectedVals[2,], apply(estimatedVals[2,,,1:2], 1, mean), 
         col=plotcolors[2,])
  
  ws.q<-apply(estimatedVals[2,,,1:2], 1, quantile, probs = c(0.025, 0.975))
  
  segments(expectedVals[2,], ws.q[1,], expectedVals[2,], ws.q[2,], 
           col=plotcolors[2,])
  
  abline(0,1)
}

###########################
#simulate community matrix#
###########################

#For debugging, to quickly generate a simulated community

# notus = 1000
# nreads = 10000
# nreps = 10
# precision = 1
# effectsize = 1.5
#  num_da = 0.25
# modeltype = "pareto_shape_4"

#Note to readers: this function is worth careful consideration. 
simCom <- function(modeltype, 
                   notus, 
                   nreads, 
                   nreps, 
                   precision, 
                   num_da = 0.25, 
                   effectsize, 
                   seed = 666,
                   thresh = 3){
  
  #Parameter explanations:
  
  #notus: number of features in composition to be modeled. Should be an even number.
  
  #nreads: the total number of observations per sample (e.g. number of reads per sample)
  
  #nreps: the number of replicates
  
  #precision: Dirichlet precision parameter. This controls how variable the draws 
  #from the Dirichlet distribution are. A higher number means more precision, and less variation 
  #among draws.
  
  #num_da: the number of taxa with different relative abundances between sampling groups. 
  #Should be an even number.
  
  #effectsize: the size of the difference in relative abundance between groups for the 
  #microbes that vary. Should be an even number.
  
  #Obtain the same seed from the input. This is to ensure we get the same 
  #deviates for each simulation replicate. 
  set.seed(seed)
  
  ####################################### 
  # Simulate Dirichlet alpha parameters #
  ####################################### 
  
  #If model type is uniform, then build multinomial sampling parameters thus:
  if(modeltype == "uniform"){
    alphas <- rep(1/notus, notus)
    #extract indices of features to apply effect size to.
    
    medium_indices <- sample(seq(1,length(alphas), by = 1),
                             size = round(length(alphas)*num_da))
    abund_indices <- NA
    rare_indices <- NA
    abund_category <- rep("med", notus)
    
    toAffect <- na.omit(c(abund_indices, medium_indices, rare_indices))
    
  }else{
    
    #This while loop is to rerun the Pareto sampling until we get a positive number of abundant taxa.
    enough <- "F"
    while(enough == "F"){
      if(modeltype == "pareto_shape_4"){
        alphas <- rpareto(notus,
                          location = 1, #floor of distribution, i.e. min. of range 
                          shape = 4)   #Degree of skew in distribution
        
        #View largest alphas if desired
        #rev(sort(alphas))[1:10]
        
        #Find indices of this vector corresponding to parameters for abundant, medium, and rare features.
        #Considered using percentiles and that didn't really seem to make sense because it does not correspond to the same thing across vectors.
        #Considering using a constant threshold (e.g. 0.1% of the data) but that didn't work when the denominator got really large, 
        #because vanishingly few features exceeded the medium and abundant thresholds.
        #So I have settled on the raw value of alpha, which is the same as a ratio between that value and the floor of the
        #Pareto distributions, because the floor was 1. 
        
        abundant <- which(alphas >= 5 ) 
        medium <- which(alphas < 5 & alphas > 2)
        rare <- which( alphas <= 2) 
        
      }else if(modeltype == "pareto_shape_0.7"){
        alphas <- rpareto(notus, 
                          location = 1, 
                          shape = 0.7)  
        
        abundant <- which(alphas >= 1000 ) 
        medium <- which(alphas < 1000 & alphas > 100)
        rare <- which( alphas <= 100) 
      }
      
      #This is to enforce a resample if we didn't get any abundant features, which happens occassionaly.
      if(length(abundant) > 0 & length(medium) > 0 & length(rare) > 0){
        
        enough <-  "T"
        
        #Create a vector describing which indices where in which abundance category. 
        #This is used in the diffAbund script to see how well the model worked in each 
        #category. 
        abund_category <- rep("NA", length(alphas))
        abund_category[abundant] <- "abund"
        abund_category[medium] <- "med"
        abund_category[rare] <- "rare"
        
        #Choose members of each abundance category for effect size application
        #NUM_DA is a proportion
        
        abund_indices <- sample(abundant, size = round(length(abundant)*(num_da / 3)))
        
        while(length(abund_indices) == 0){
          abund_indices <- sample(abundant, size = 1)
        }
        
        medium_indices <- sample(medium, size = round(length(medium)*(num_da / 3)))
        
        while(length(medium_indices) == 0){
          medium_indices <- sample(medium, size = 1)
        }
        
        rare_indices <- sample(rare, size = round(length(rare)*(num_da / 3)))
        
        while(length(rare_indices) == 0){
          rare_indices <- sample(rare, size = 1)
        }
        
        toAffect <- na.omit(c(abund_indices, medium_indices, rare_indices))
        
        #Adding this to avoid odd ball situations where toAffect is very small. 
        if(length(toAffect) < 10){
          enough <- "F"
        }
      }
    }
  }
  
  #Need to make sure the exact same features differ from one group to the next and the sum of each vector is the same.
  #To do this, we sample indices from each of the feature abundance classes (abund, med, rare)
  #and duplicate those features by appending them onto the end of the alpha vector.
  #Then, I multiply those appended features by the effect size, this makes alphas for group 1, 
  #Next, I multiply the features in group 2 that were not appended by the effect size.
  #Thus the last features in the alpha vector differ among groups and the features corresponding to the 
  #indices in toAffect also differ.
  
  #For example: 
  #Original alpha vector: 0.5,0.1,0.01
  #New alpha vector: 0.5,0.1,0.01,  0.5,0.1,0.01
  
  #After applying effect sizes and making vectors for both groups: 
  #I added tabs to make this a bit clearer to read.
  
  #group 1: 0.5,0.1,0.01,       0.5*effect,0.1*effect,0.01*effect
  #group 2: 0.5*effect,0.1*effect,0.01*effect,      0.5,0.1,0.01
  
  #Apply effect size to make group 1 alpha parameter vector
  alphas_group1 <- c(alphas, alphas[toAffect]*effectsize)
  
  #Apply effect size to make group 2 alpha parameter vector
  alphas_group2 <- c(alphas, alphas[toAffect])
  alphas_group2[toAffect] <- alphas_group2[toAffect]*effectsize
  
  abund_category <- c(abund_category, abund_category[toAffect])
  
  #Sanity checks 
  #alphas_group2 - alphas_group1
  #Make sure vectors add up to the same thing
  sum(alphas_group2)
  sum(alphas_group1)
  
  #Define a vector of indices for features that should differ between groups. Note 
  #That this is for labeling purposes and gets overwritten after sampling from 
  #Dirichlet and multinomial distributions. 
  
  differentTaxa <- c(toAffect, 
                     seq(length(alphas) + 1,
                         length(alphas) + length(toAffect), by = 1) 
  )
  
  #Sample from the Dirichlet to make p parameters for the multinomial. 
  #This is the part that should vary among replicates, therefore we redefine the seed. 
  
  set.seed(NULL)
  
  comMat <- matrix(0, 
                   ncol = length(alphas_group1), #Both groups are same dimensions
                   nrow = nreps)
  
  for(i in 1:(nreps / 2)){
    comMat[i,] <- rmultinom(1, 
                            nreads, 
                            prob = rdirichlet(1, alpha = alphas_group1*precision) )
    
  }
  for(i in (1 + (nreps / 2)):nreps){
    comMat[i,] <- rmultinom(1, 
                            nreads, 
                            prob = rdirichlet(1, alpha = alphas_group2*precision))
    
  }
  
  #Assign meaningful names to the different features,
  #This will allow us to identify which features should differ later on. 
  colnames(comMat) <- paste("Feature_", 
                            seq(1, length(alphas_group1), 
                                by = 1), 
                            abund_category, sep = "")
  
  colnames(comMat)[differentTaxa] <- paste(colnames(comMat)[differentTaxa], 
                                           "_different", 
                                           sep = "")
  comMat <- data.frame(comMat)
  
  #Remove features with zero counts / low abundance across replicates.
  #I choose three reads as the default somewhat arbitrarily. However, it is common for singletons, doubletons to be removed
  #from sequencing data as they could be spurious (i.e. PCR errors, etc.). Therefore, three seems
  #a practical threshold. 
  
  #Remove these features from the differentTaxa object as well by remaking it.
  #This is inelegant, but gets the job done.
  
  #comMat <- comMat[, which(colSums(comMat) > thresh)]
  differentTaxa <- grep("different", colnames(comMat))
  #abund_category <- abund_category[which(colSums(comMat) > thresh )]
  
  #Save alphas as proportions. We use this to compute model RMSE
  #NOTE: we do not include low abundance features here, see following lines.
  #comprop <- rbind(alphas_group1[which(colSums(comMat) > thresh )] / sum(alphas_group1[which(colSums(comMat) > thresh )]),
  #               alphas_group2[which(colSums(comMat) > thresh )] / sum(alphas_group2[which(colSums(comMat) > thresh )]))
  
  #Remove uncommon features
  #Remove those with < 100 reads
  #Remove those that don't occur in 15% or more of samples
  
  #Add a 1 to avoid infinite density slicer error generation by rjags
  #This can still be an issue so long as there is a zero anywhere in the simulated data.
  comMat <- comMat + 1
  
  return(list(simulatedCommunity = comMat, 
              differentTaxa = differentTaxa, 
              Parameters_alpha1_2_pi1_2 = rbind(alphas_group1/sum(alphas_group1),
                                                alphas_group2/sum(alphas_group2)),
              abund_category = abund_category,
              exp_abund_differ = grep("abund_different", colnames(comMat)),
              exp_med_differ = grep("med_different", colnames(comMat)),
              exp_rare_differ = grep("rare_different", colnames(comMat))
  ))
}

###########################################################################
#splitter - to split a matrix by row and turn each row into a list element#
###########################################################################

splitter <- function(x){split(x, seq(nrow(x)))}

######################################################
#Wilcoxon tests: conduct Mann-Whitney U/Wilcoxon test#
######################################################

#this loop runs a test for each taxon/element in the compositions
wilcoxonTester <- function(simulatedData = simcom_out$simulatedCommunity,
                           groupings, 
                           otus, 
                           fname,
                           abund_category = simcom_out$abund_category,
                           Parameters_alpha1_2_pi1_2 = simcom_out$Parameters_alpha1_2_pi1_2,
                           different = simcom_out$differentTaxa,
                           rare_different = simcom_out$exp_rare_differ,
                           med_different = simcom_out$exp_med_differ,
                           abund_different = simcom_out$exp_abund_differ){
  
  testOut <- list()
  
  group1start <- min(which(groupings == unique(groupings)[1]))
  group1end <-  max(which(groupings == unique(groupings)[1]))
  group2start <- min(which(groupings == unique(groupings)[2]))
  group2end <-  max(which(groupings == unique(groupings)[2]))
  
  for(i in 1:dim(simulatedData)[2]){
    testOut[[i]] <- wilcox.test(simulatedData[group1start:group1end,i], 
                                simulatedData[group2start:group2end,i])
  }
  
  positives <- which(lapply(testOut,FUN=function(x){x$p.value <= 0.05})==TRUE)
  negatives <- which(lapply(testOut,FUN=function(x){x$p.value > 0.05})==TRUE)
  
  positivesBH <- which(p.adjust(unlist(lapply(testOut, function(x){x[c('p.value')]})),"fdr")<=0.05)
  negativesBH <- which(p.adjust(unlist(lapply(testOut, function(x){x[c('p.value')]})),"fdr")>0.05)
  
  positives_diff <- data.frame(positives, abund_category[positives])
  negatives_diff <- data.frame(negatives, abund_category[negatives])
  
  positives_diffBH <- data.frame(positivesBH, abund_category[positivesBH])
  negatives_diffBH <- data.frame(negativesBH, abund_category[negativesBH])
  
  # fp_effects <- Parameters_alpha1_2_pi1_2[1, positives[!(positives %in% different)]] - Parameters_alpha1_2_pi1_2[2, positives[!(positives %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/", fname, "_fp_effects_wilcox.csv", sep = ""))
  # 
  # fp_effects <- Parameters_alpha1_2_pi1_2[1, positivesBH[!(positivesBH %in% different)]] - Parameters_alpha1_2_pi1_2[2, positivesBH[!(positivesBH %in% different)]]
  # 
  # write.csv(fp_effects, file = paste("/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/output/", fname, "_fp_effects_wilcoxBH.csv", sep = ""))
  # 
  return(
    list(
      performanceAll = fpr(
        otus = otus,
        xpos = positives,
        xneg = negatives,
        different = different
      ),
      performanceBH = fpr(
        otus = otus,
        xpos = positivesBH,
        xneg = negativesBH,
        different = different
      ),
      performance_abund = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "abund")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "abund")],
        different = abund_different 
      ),
      performance_med = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "med")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "med")],
        different = med_different
      ),
      performance_rare = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. == "rare")],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. == "rare")],
        different = rare_different
      ),
      performance_abundBH = fpr(
        otus = otus[which(abund_category == "abund")],
        xpos = positives_diffBH$positives[which(positives_diffBH$abund_category.positivesBH. == "abund")],
        xneg = negatives_diffBH$negatives[which(negatives_diffBH$abund_category.negativesBH. == "abund")],
        different = abund_different
      ),
      performance_medBH = fpr(
        otus = otus[which(abund_category == "med")],
        xpos = positives_diffBH$positives[which(positives_diffBH$abund_category.positivesBH. == "med")],
        xneg = negatives_diffBH$negatives[which(negatives_diffBH$abund_category.negativesBH. == "med")],
        different = med_different
      ),
      performance_rareBH = fpr(
        otus = otus[which(abund_category == "rare")],
        xpos = positives_diffBH$positives[which(positives_diffBH$abund_category.positivesBH. == "rare")],
        xneg = negatives_diffBH$negatives[which(negatives_diffBH$abund_category.negativesBH. == "rare")],
        different = rare_different
      ),
      performance_med_hi = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diff$positives[which(positives_diff$abund_category.positives. %in% c("abund", "med"))],
        xneg = negatives_diff$negatives[which(negatives_diff$abund_category.negatives. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      ),
      performance_med_hiBH = fpr(
        otus = otus[which(abund_category %in% c("abund", "med"))],
        xpos = positives_diffBH$positives[which(positives_diffBH$abund_category.positivesBH. %in% c("abund", "med"))],
        xneg = negatives_diffBH$negatives[which(negatives_diffBH$abund_category.negativesBH. %in% c("abund", "med"))],
        different = c(med_different, abund_different)
      )
    )
  )
}


#################################
#Conduct simulation experiment ##
#################################

##############################################################
#Run the main function after pulling parameter line from STDIN#
##############################################################

#Recall that in the wrap*pl script different numbers are substituted into STDIN. 
#This lets us read different lines of the input file via this function.
runner <- function() {
  inargs <- commandArgs(trailingOnly = TRUE) 
  print("Iteration number:")
  print(inargs)
  #Cuts out the calls to R that Rscript includes automatically, eg --slave --no-restore, etc
  parameters <- read.csv(#"./serialRunParams.csv")
    "/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/serialRunParamsv2.csv", stringsAsFactors = F)
  
  tryCatch(main(as.vector(unlist(parameters[inargs[1],]))),
           warnings=function(w){
             fileConn<-file(paste(#"./warnings/",
               "/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/warnings/warnings",
               parameters[1],
               "otus",
               parameters[2],
               "reads",
               parameters[3],
               "reps",
               parameters[4],
               "theta",
               parameters[5],
               "diffabund",
               parameters[6],
               "alpha",
               parameters[7],
               "effectsize",
               parameters[8],
               "simulationrep.csv",
               sep = ""))
             
             writeLines(names(warnings()), fileConn)
             on.exit(close(fileConn))
             return(NA)
           },
           errors=function(e){
             fileConn<-file(paste(#"./errors/",
               "/project/microbiome/users/jharri62/DirichletSimulationPaper/dmmodelvetting/errors/errors",
               parameters[1],
               "otus",
               parameters[2],
               "reads",
               parameters[3],
               "reps",
               parameters[4],
               "theta",
               parameters[5],
               "diffabund",
               parameters[6],
               "alpha",
               parameters[7],
               "effectsize",
               parameters[8],
               "simulationrep.csv",
               sep = ""))
             
             writeLines(names(errors()), fileConn)
             on.exit(close(fileConn))
             return(NA)
           }
  )
}

runner()

# # # #debugging code
parameters <- read.csv("./redo.csv", stringsAsFactors = F)

parameters <- as.vector(unlist(parameters[2,]))
#  main(parameters)
#parameters <- subset(parameters, X1 == "pareto_shape_0.7")

