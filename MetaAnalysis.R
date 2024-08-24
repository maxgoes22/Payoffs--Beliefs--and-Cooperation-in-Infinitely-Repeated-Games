####################
#
# Main
#
#######################


#Working Directory
wd <- "C:/Users/Maximilian Andres/Nextcloud2/TextMiningProject/05_Gain_and_Loss/Data Analysis/"
wd <- "/Users/juri/Nextcloud/TextMiningProject/05_Gain_and_Loss/Data Analysis/"

#Set Directory of the Meta Data
filename <- "Meta Data/Data_DalBoFrechette2018.txt"

#source functions
source(paste(wd,"MetaFunctions.R",sep=""))

#load data
metadata <- load_data(wd, filename)

# number of observations
nrow(metadata)
nrow(metadata[metadata$round >1, ])
nrow(metadata[metadata$round == 1, ])

# summary statistics
round(min(metadata[metadata$round == 1, ]$g), 2)
round(max(metadata[metadata$round == 1, ]$g), 2)
round(min(metadata[metadata$round == 1, ]$l), 2)
round(max(metadata[metadata$round == 1, ]$l), 2)
round(cor(metadata[metadata$round == 1, ]$g, metadata[metadata$round == 1, ]$l), 2)

# belief construction
ranges <- c(-1,0.2,0.4,0.6,0.8,1)
lambda <- c(0.3, 0.4, 0.6, 0.7, 0.5)

#install.packages("doParallel")
library("doParallel")
cl = makeCluster(ceiling(detectCores()*0.7))
registerDoParallel(cl)
for (i in 1:length(lambda)){
  metadata[paste("Markov", lambda[i], sep="")] <- belief_construction_markov_parallel(metadata, lambda = lambda[i], varofinterest = "CooperationRate")
  print(paste("Done with run ", i, " of ", length(lambda), sep = ""))
}
stopCluster(cl)


# plot histogram of the constructed beliefs
plot_histogram(metadata, lambda = 0.5, appendix=TRUE)
ggsave(paste(wd,"Figures Meta Data/Figure_Histogram_Markov",0.5,"_Main",".pdf",sep=''), width = 3.5, height = 3.5)

for (lambdavariable in lambda) {
  plot_histogram(metadata, lambda = lambdavariable, appendix=FALSE)
  ggsave(paste(wd,"Figures Meta Data/Figure_Histogram_Markov",lambdavariable,"_Appendix",".pdf",sep=''), width = 3.5/1.5, height = 3.5/1.5)
}

# plot confusion matrix and print the accuracy of the constructed beliefs
plot_confusion_matrix(metadata, lambda = 0.5, appendix=FALSE)
ggsave(paste(wd,"Figures Meta Data/Figure_Confusion_Markov",0.5,"_Main",".pdf",sep=''), width = 3.5, height = 3.5)

for (lambdavariable in lambda) {
  plot_confusion_matrix(metadata, lambda = lambdavariable, appendix=TRUE)
  ggsave(paste(wd,"Figures Meta Data/Figure_Confusion_Markov",lambdavariable,"_Appendix",".pdf",sep=''), width = 3.5/1.5, height = 3.5/1.5)
}

# the accuracy of the sub-game perfect critical discount factor
round(nrow(metadata[metadata$round == 1 & metadata$sgpeDummy == metadata$Cooperation, ]) / nrow(metadata[metadata$round == 1, ]), 2)

# the accuracy of the risk-dominant critical discount factor
round(nrow(metadata[metadata$round == 1 & metadata$rdDummy == metadata$Cooperation, ]) / nrow(metadata[metadata$round == 1, ]), 2)

# Partial Effects 
MarkovRegressionData <- partial_effect_regression(metadata = metadata[metadata$round == 1,], lambda = lambda, ranges = ranges, beliefofinterest = "Markov")
table_regression_data(regressiondata = MarkovRegressionData, lambda = 0.5, asLaTeX = F)

# Gain
plot_layers(MarkovRegressionData[MarkovRegressionData$LearningParameter == 0.5,], ranges, layer = "g", colour_learning = c("#d95f02"), y_axis_title ="Partial Effect of the Gain on Cooperation", yrange = 0.25)
ggsave(paste(wd,"Figures Meta Data/Figure_PartialEffects_Gain_Markov05.pdf",sep=""), width = 2.5, height = 3.5)

# Loss
plot_layers(MarkovRegressionData[MarkovRegressionData$LearningParameter == 0.5,], ranges, layer = "l", colour_learning = c("#7570b3"), y_axis_title ="Partial Effect of the Loss on Cooperation", yrange = 0.25)
ggsave(paste(wd,"Figures Meta Data/Figure_PartialEffects_Loss_Markov05.pdf",sep=""), width = 2.5, height = 3.5)

# Difference
plot_layers(MarkovRegressionData[MarkovRegressionData$LearningParameter == 0.5,], ranges, layer = "diff", colour_learning = c("#1b9e77"), y_axis_title ="Effect of the Gain - Effect of the Loss", yrange = 0.35)
ggsave(paste(wd,"Figures Meta Data/Figure_PartialEffects_Diff_Markov05.pdf",sep=""), width = 2.5, height = 3.5)

# Gain All
plot_layers(MarkovRegressionData, ranges, layer = "g", colour_learning = c("#d95f02", "#d95f02", "#d95f02", "#d95f02", "#d95f02"), y_axis_title ="Partial Effect of the Gain on Cooperation", yrange = 0.25)
ggsave(paste(wd,"Figures Meta Data/Figure_PartialEffects_Gain_MarkovAll.pdf",sep=""), width = 2.5, height = 3.5)

# Loss All
plot_layers(MarkovRegressionData, ranges, layer = "l", colour_learning = c("#7570b3", "#7570b3", "#7570b3", "#7570b3", "#7570b3"), y_axis_title ="Partial Effect of the Loss on Cooperation", yrange = 0.25)
ggsave(paste(wd,"Figures Meta Data/Figure_PartialEffects_Loss_MarkovAll.pdf",sep=""), width = 2.5, height = 3.5)

# Difference All
plot_layers(MarkovRegressionData, ranges, layer = "diff", colour_learning = c("#1b9e77", "#1b9e77", "#1b9e77", "#1b9e77", "#1b9e77"), y_axis_title ="Effect of the Gain - Effect of the Loss", yrange = 0.35)
ggsave(paste(wd,"Figures Meta Data/Figure_PartialEffects_Diff_MarkovAll.pdf",sep=""), width = 2.5, height = 3.5)

#plot learning path
for (lambdavariable in lambda) {
  plot_learning_path(lambdavariable)
  ggsave(paste(wd,"Figures Meta Data/Figure_SimMarkov_PofC_lambda_of_", lambdavariable,".pdf", sep = ""), width = 2.5, height = 3.5)
}

