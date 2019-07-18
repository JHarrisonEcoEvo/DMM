# DMM

The scripts provided in this repo allow reproduction of the simulation experiment conducted in the manuscript, "Dirichlet-multinomial modelling outperforms alternatives for analysis of microbiome and other ecological count data" by Joshua G. Harrison, W. John Calder, Vivaswat Shastry, and C. Alex Buerkle. 

For the impatient: the SimulationCode.R script is what we used to conduct our simulation experiment. 

Description of files (alphabetical): 

competeModelPlot.R - script to create panels in Figures 4 & S3.

data/allv4.csv – output from simulation experiment. File created by SimulationCode.R. This is wrangled in plotting scripts.

DM.stan – Stan DMM model

effectSizeBoxplot.R – script to creat Figures 3, S1, S2, and panel b in Fig. 5.

lm_of_simulation_params.R – script to perform linear modeling analysis of the influence of data attributes on model performance. Results are shown in Tables S2 & S3.

mcmcComparisonPlot.R – script to make Figure 5.

README.md – a highly entertaining document!

reanalysisLungsHMC.R – analysis of a portion of the data published in Duvallet et al 2019. See main text. This script requires "patient_clinical_metadata.csv" and the directory "rosen_mincount10_maxee2_trim200_results_forpaper" from Duvallet et al. These data can be accessed at: https://doi.org/10.5281/zenodo.2678107.
If you use these data please cite Duvallet et al. 2019.

simParamMaker.R – a simple script to make a file that is read by SimulationCode.R

SimulationCode.R – the most important file here. This is the code used to conduct the simulation experiment. It is a collection of functions that are presented alphabetically. Each function is called from a "main" function, which is the place to start for those perusing this script. 