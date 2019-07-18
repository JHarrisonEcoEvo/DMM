
library(gtools)
library(rstan)
library(rjags)
library(shinystan)
library(VGAM) #Pareto
set.seed(1245)
options(scipen=99)


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

#Function to calculate initial values for MCMC chains
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
      pi[i,,j] <- (apply(dat[starts[i]:ends[i],1:notusinit], 2, mean) + 1) /
        (sum(apply(dat[starts[i]:ends[i],1:notusinit], 2, mean) + 1))
      pi[i,,j] <- pi[i,,j]
    }
  }
  
  #Generate initial values for theta
  theta <- matrix(nrow = 2,
                  ncol = length(starts))
  for(j in 1:2){
    for(i in 1:length(starts)){
      thetas <- alphas / pi[i,,1] #There is no reason to use both parts of array
      
      theta <- mean(thetas[which(thetas < 15000)])
      
      #Here we avoid initializing using high values of theta that seem 
      #unrealistic.
    }
  }
  
  return(list(list("p" = p[,,1],
                   "pi" = pi[,,1]
                  # "alpha" =  c(alphas)/sum(alphas),
                   #"theta"  = theta
                  ),
              list("p" = p[,,2],
                   "pi" = pi[,,2]
                   #"alpha" =  c(alphas2)/sum(alphas2),
                   #"theta"  = theta
                   )
  )
  )
}


# Function to simulate data.
simCom <- function(modeltype,
                   notus,
                   nreads,
                   nreps,
                   precision,
                   num_da = 0.25,
                   effectsize,
                   seed = 666,
                   thresh = 3){
  set.seed(seed)
  
  #######################################
  # Simulate Dirichlet alpha parameters #
  #######################################
  enough <-  "F"
  while(enough == "F"){
    alphas <- VGAM::rpareto(notus,
                            scale = 1, #floor of distribution, i.e. min. of range
                            shape = 4)   #Degree of skew in distribution
    
    abundant <- which(alphas >= 5 )
    medium <- which(alphas < 5 & alphas > 2)
    rare <- which( alphas <= 2)
    
    #This is to enforce a resample if we didn't get any abundant features
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
  
  #Must make sure the exact same features differ from one group to the next and
  #the sum of each vector is the same.
  #To do this, we sample indices from each of the feature abundance classes
  #(abund, med, rare) and duplicate those features by appending them onto the
  #end of the alpha vector.
  
  #Then, I multiply those appended features by the effect size, this makes alphas for
  #group 1,
  #Next, I multiply the features in group 2 that were not appended by the effect size.
  #Thus the last features in the alpha vector differ among groups and the features
  #corresponding to the indices in toAffect also differ.
  
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
  
  #Remove these features from the differentTaxa object as well by remaking it.
  #This is inelegant, but gets the job done.
  
  comMat <- comMat[, which(colSums(comMat) > thresh )]
  differentTaxa <- grep("different", colnames(comMat))
  
  #Save alphas as proportions. We use this to compute model RMSE
  #NOTE: we do not include low abundance features here, see following lines.
  comprop <- rbind(alphas_group1[which(colSums(comMat) > thresh )] /
                     sum(alphas_group1[which(colSums(comMat) > thresh )]),
                   alphas_group2[which(colSums(comMat) > thresh )] /
                     sum(alphas_group2[which(colSums(comMat) > thresh )]))
  
  #Add a 1 to avoid infinite density slicer error generation by rjags
  comMat <- comMat + 1
  
  return(list(simulatedCommunity = comMat,
              differentTaxa = differentTaxa,
              Parameters_alpha1_2_pi1_2 = comprop,
              abund_category = abund_category,
              exp_abund_differ = grep("abund_different", colnames(comMat)),
              exp_med_differ = grep("med_different", colnames(comMat)),
              exp_rare_differ = grep("rare_different", colnames(comMat))
  ))
}

simcom_out <- simCom(modeltype = "pareto_shape_4", 
                     notus = 100, #Note that if this is too small, this function will hang bc it is not generating enough features in each abundance class
                     nreads = 1000, 
                     nreps = 40, 
                     precision = 1, 
                     num_da = 0.25, 
                     effectsize = 2, 
                     seed = 666,
                     thresh = 3)

#Make a design variable, describing which sample goes with which group. 

groupings <- c(rep("group1", 40 / 2),
               rep("group2", 40 / 2))

DM<-stan_model("DM.stan", model_name="DM")

x <- proc.time()

inits <- initcalculator(dat = simcom_out$simulatedCommunity,
                        starts = c(min(which(groupings == unique(groupings)[1])),
                                   min(which(groupings ==unique(groupings)[2]))),
                        ends = c(max(which(groupings == unique(groupings)[1])),
                                 max(which(groupings == unique(groupings)[2]))))
fitstan_VB <-vb(DM, 
                data=list("datamatrix"=simcom_out$simulatedCommunity, 
                          "nreps"=nrow(simcom_out$simulatedCommunity), 
                          "notus"=ncol(simcom_out$simulatedCommunity),
                          "N"=2, 
                          "start" = c(min(which(groupings == unique(groupings)[1])),
                                      min(which(groupings ==unique(groupings)[2]))),
                          "end" = c(max(which(groupings == unique(groupings)[1])),
                                    max(which(groupings == unique(groupings)[2])))), 
                algorithm="meanfield", 
                output_samples=500, 
                check_data = T,
                seed=123,
                pars<-c("pi")
                
          
                ,init = list("p"= list(inits[[1]]$p, inits[[2]]$p),
                             "pi"= inits[[1]]$pi)
) 

print("Time for Stan to run using VB algorithm")
vb_time <- proc.time()-x

#Run NUTS algorithm
x <- proc.time()
fitstan_NUTS<-sampling(DM, 
                       data=list(datamatrix=simcom_out$simulatedCommunity, 
                                 nreps=nrow(simcom_out$simulatedCommunity), 
                                 notus=ncol(simcom_out$simulatedCommunity), N=2, 
                                 start = c(min(which(groupings == unique(groupings)[1])),
                                           min(which(groupings ==unique(groupings)[2]))),
                                 end = c(max(which(groupings == unique(groupings)[1])),
                                         max(which(groupings == unique(groupings)[2])))), 
                       chains=2, 
                       control = list(max_treedepth = 15),
                       warmup=10, 
                       iter=20, 
                       thin=2, 
                       algorithm="NUTS", 
                       cores=2, 
                       pars<-c("pi"),
                       init = initcalculator(dat = simcom_out$simulatedCommunity, 
                                             starts = c(min(which(groupings == unique(groupings)[1])),
                                                        min(which(groupings == unique(groupings)[2]))), 
                                             ends = c(max(which(groupings == unique(groupings)[1])),
                                                      max(which(groupings == unique(groupings)[2])))),
                       verbose=T) 
print("Time for Stan to run using NUTS algorithm")
nutsTime <- proc.time()-x



community.model.level <- "model{
for(i in 1:N){
for(j in start[i]:end[i]){
datamatrix[j,] ~ dmulti(p[j,], nreads[j])
p[j,1:notus] ~ ddirch(pi[i,]*theta[i])
}

pi[i,1:notus] ~ ddirch(alpha)
theta[i] ~ dunif(0.1, 4000)
}

for(k in 1:notus){
alpha[k] <- 0.0000001
}
}"

#########################
# Function to run model #
#########################

modelRun <- function(SimulatedData, 
                     groupings){
  
  #Specify model
  sim.mod.jags <- rjags::jags.model(textConnection(community.model.level),
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
  
  #Adapt model for as long as needed
  iter_needed <- 0
  y = FALSE
  while(y == FALSE){
    y <-  adapt(sim.mod.jags, 
                n.iter = 10,
                end.adaptation = FALSE)
    iter_needed <- 10 + iter_needed
    if(iter_needed > 50){break} 
  }
  
  #Burn in model. This should lead to convergence on the high density interval 
  #of posterior parameter space. 
  
  update(sim.mod.jags, 
         n.iter = 50) 
  
  #Here we extract the MCMC samples 
  sim.mod.sam <- rjags::jags.samples(model = sim.mod.jags, 
                                     variable.names=c("pi"
                                                      ),
                                     n.iter = 40, thin = 4)
  return(sim.mod.sam)
}

sim.mod.sam <- modelRun(SimulatedData = simcom_out$simulatedCommunity, 
                        groupings = groupings) 

est.pi <- extract(fitstan_NUTS,"pi")
est.pi_vb <- extract(fitstan_VB,"pi")

#Plot

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#order by expected relative abundance
orderUse <- order(simcom_out$Parameters_alpha1_2_pi1_2[1,])

par(oma = c(1,1,0,0))
plot(apply(est.pi$pi[,1,], 2, FUN=mean)[orderUse], 
     pch = 16,
     ylab = "",
     xlab = "Feature (ordered by relative abundance)",
     frame.plot = F,
     yaxt = "n",
     xaxt = "n",
     ylim = c(0.005, .045),
     xlim = c(0,110),
     col = add.alpha("darkslateblue", 0.6))

axis(side = 2,at = c(seq(0.005, 0.045, by = 0.01)), labels = c(seq(0.005, 0.045, by = 0.01)),
     las = 2)
axis(side = 1, 
     at = c(seq(0, 100, by = 20),110), 
     labels = c(seq(0, 100, by = 20),110))

text(expression(paste(pi, " estimate"), sep = ""),
     x = -19, y = 0.025, srt = 90,xpd = NA, cex = 1.5)

points(fitstan_VB@sim$est$pi[1,orderUse],
       col = add.alpha("darkseagreen4", 0.6),
       pch = 16)

points(apply(sim.mod.sam$pi[1,orderUse,,1:2], 1, mean),
       col = add.alpha("coral3", 0.6),
       pch = 16)

points(simcom_out$Parameters_alpha1_2_pi1_2[1,orderUse],
       col = "black",
       pch = 3,
       cex = 0.3)

#jags
segments(x0 = c(seq(1, 110, by = 1) - 0.3),
         x1 = c(seq(1, 110, by = 1) + 0.3),
         y0 = apply(sim.mod.sam$pi[1,orderUse,,1:2], 1, HDIofMCMC, credMass = 0.95)[1,],
         y1 = apply(sim.mod.sam$pi[1,orderUse,,1:2], 1, HDIofMCMC, credMass = 0.95)[1,],
       col = add.alpha("coral3", 0.6))

segments(x0 = c(seq(1, 110, by = 1) - 0.3),
         x1 = c(seq(1, 110, by = 1) + 0.3),
         y0 = apply(sim.mod.sam$pi[1,orderUse,,1:2], 1, HDIofMCMC, credMass = 0.95)[2,],
         y1 = apply(sim.mod.sam$pi[1,orderUse,,1:2], 1, HDIofMCMC, credMass = 0.95)[2,],
         col = add.alpha("coral3", 0.6))

segments(x0 = c(seq(1, 110, by = 1) - 0.05),
         x1 = c(seq(1, 110, by = 1) - 0.05),
         y0 = apply(sim.mod.sam$pi[1,orderUse,,1:2], 1, HDIofMCMC, credMass = 0.95)[1,],
         y1 = apply(sim.mod.sam$pi[1,orderUse,,1:2], 1, HDIofMCMC, credMass = 0.95)[2,],
         col = add.alpha("coral3", 0.6),
         lty = 3)

#nuts
segments(x0 = c(seq(1, 110, by = 1) - 0.3),
         x1 = c(seq(1, 110, by = 1) + 0.3),
         y0 = apply(est.pi$pi[,1,orderUse], 2, FUN=HDIofMCMC)[1,],
         y1 = apply(est.pi$pi[,1,orderUse], 2, FUN=HDIofMCMC)[1,],
         col = add.alpha("darkseagreen4", 0.6),
         pch = 1)

segments(x0 = c(seq(1, 110, by = 1) - 0.3),
         x1 = c(seq(1, 110, by = 1) + 0.3),
         y0 = apply(est.pi$pi[,1,orderUse], 2, FUN=HDIofMCMC)[2,],
         y1 = apply(est.pi$pi[,1,orderUse], 2, FUN=HDIofMCMC)[2,],
         col = add.alpha("darkseagreen4", 0.6),
         pch = 1)

segments(x0 = c(seq(1, 110, by = 1) + 0.05),
         x1 = c(seq(1, 110, by = 1) + 0.05),
         y0 =apply(est.pi$pi[,1,orderUse], 2, FUN=HDIofMCMC)[1,],
         y1 = apply(est.pi$pi[,1,orderUse], 2, FUN=HDIofMCMC)[2,],
         col = add.alpha("darkseagreen4", 0.6),
         lty = 3)

#vb
segments(x0 = c(seq(1, 110, by = 1) - 0.3),
         x1 = c(seq(1, 110, by = 1) + 0.3),
         y0 = apply(est.pi_vb$pi[,1,orderUse], 2, FUN=HDIofMCMC)[1,],
         y1 = apply(est.pi_vb$pi[,1,orderUse], 2, FUN=HDIofMCMC)[1,],
         col = add.alpha("darkslateblue", 0.6),
         pch = 1)

segments(x0 = c(seq(1, 110, by = 1) - 0.3),
         x1 = c(seq(1, 110, by = 1) + 0.3),
         y0 = apply(est.pi_vb$pi[,1,orderUse], 2, FUN=HDIofMCMC)[2,],
         y1 = apply(est.pi_vb$pi[,1,orderUse], 2, FUN=HDIofMCMC)[2,],
         col = add.alpha("darkslateblue", 0.6),
         pch = 1)

segments(x0 = c(seq(1, 110, by = 1)),
         x1 = c(seq(1, 110, by = 1)),
         y0 =apply(est.pi_vb$pi[,1,orderUse], 2, FUN=HDIofMCMC)[1,],
         y1 = apply(est.pi_vb$pi[,1,orderUse], 2, FUN=HDIofMCMC)[2,],
         col = add.alpha("darkslateblue", 0.6),
         lty = 3)
