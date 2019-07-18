#J. Harrison
#University of Wyoming

#specify test parameters and turn them into a list
modeltype <- c("uniform", "pareto_shape_0.7", "pareto_shape_4")
notus <- c(500,2000)
nreads <- c(50000,10000)
nreps <- c(10,50)
theta <- c(0.5,3) #precision
diffabund <- c(0.25)
effect <- c(1.1, 1.5, 2) #Base value needs to be 1. 
simreps <- 1

#expand grid 
parameters <- expand.grid(as.character(modeltype),
                          notus,
                 nreads, 
                 nreps, 
                 theta, 
                 diffabund,
                 effect
                 )

i <- 1
parametersReplicated <- parameters
while(i < simreps){
  parametersReplicated <- rbind(parametersReplicated, parameters)
  i <- i + 1
}

splitter <- function(x){split(x, seq(nrow(x)))}
parametersReplicated[,1] <- as.character(parametersReplicated[,1] )

inparam <- splitter(parametersReplicated)


params <- matrix(unlist(inparam), 
                  ncol = dim(parameters)[2], 
                  byrow = T)

#add index to reference simreps. This is to prevent overwriting output files
indices <- seq(1, length(params[,1]), by = 1)

write.csv(data.frame(params,indices) , 
           row.names = F, 
           file = "serialRunParamsv2.csv")
