#Packages
#install.packages("tidyverse")
library("tidyverse")
#install.packages("mfx")
library("mfx")
#install.packages("foreach")
library("foreach")
#install.packages("doParallel")
library("doParallel")

#Working Directory
wd <- "C:/Users/Maximilian Andres/Nextcloud2/TextMiningProject/05_Gain_and_Loss/Data Analysis/"
#wd <- "/Users/juri/Nextcloud/TextMiningProject/05_Gain_and_Loss/Data Analysis/"
filename <- "MetaStudy/Data_DalBoFrechette2018.txt"


##################################
#
# Functions
#
##################################

#load the data
load_data <- function(wd, filename){
  
  #read meta data
  metadata <- read.delim(paste(wd, filename, sep=""), header = TRUE, sep = "\t", dec = ".")
  
  #just focus on infinitely repeated games
  metadata <- metadata[metadata$delta >0,]

  #metrices
  metadata$sizeBAD <- ((1-metadata$delta)*(metadata$p-metadata$s))/((metadata$r-metadata$p)-(1-metadata$delta)*(metadata$t+metadata$s-(2*(metadata$p))))
  metadata$sgpe <- ((metadata$t-metadata$r)/(metadata$t-metadata$p))
  metadata$rd <- ((metadata$t-metadata$r + metadata$p-metadata$s)/(metadata$t-metadata$s))
  
  # is grim subgame perfect?
  metadata$sgpeDummy = ifelse(metadata$delta >= metadata$sgpe, 1, 0)
  
  # is grim risk dominant?
  metadata$rdDummy = ifelse(metadata$delta >= metadata$rd, 1, 0)
  
  #cooperation variable binary
  metadata$Cooperation  <- ifelse(metadata$coop == "coop",1,0)
  
  #other cooperation variable binary
  metadata$OtherCooperation  <- ifelse(metadata$ocoop == "coop",1,0)
  
  #mutual cooperation variable binary
  metadata$MutualCooperation  <- ifelse(metadata$coop == "coop" & metadata$ocoop == "coop",1,0)
  
  #cooperation rate variable
  metadata$CooperationRate  <- (metadata$Cooperation + metadata$OtherCooperation)/(2)
  
  #payoff variable
  metadata$Payoff  <- ifelse((metadata$Cooperation == 1) & (metadata$OtherCooperation == 1),1,
                              ifelse((metadata$Cooperation == 1) & (metadata$OtherCooperation == 0),-metadata$l,
                                     ifelse((metadata$Cooperation == 0) & (metadata$OtherCooperation == 1),1+metadata$g,
                                            ifelse((metadata$Cooperation == 1) & (metadata$OtherCooperation == 0),1,0))))
  
  #variable for partial effect regressions
  metadata$gplusl <- metadata$g +metadata$l
  
  #in the meta data set of Dal Bo and Frechette (2018), for subjects in the Duffy and Ochs (2009) paper, the subject id's are not unique. Thus, we drop the data of this paper for our analysis.
  metadata <- metadata[metadata$paper != "Duffy and Ochs 2009",]
  
  return(metadata)
  
}

#construct belief Markov Chain considering the rounds (PARALLEL)
belief_construction_markov_parallel_B1 <- function(metadata, lambda, varofinterest){
  
  metadata$Output <- NA
  prior <- 0.5
  id <- unique(metadata$id)
  
  # Parallel computation
  result <- foreach(i = 1:length(id), .combine = rbind, .packages = c("tidyverse")) %dopar% {
    
    localdata <- metadata[metadata$id == id[i], ]
    
    supergame <- sort(unique(localdata$supergame))
    
    for (s in seq_along(supergame)) {
      
      period <- sort(unique(localdata[localdata$supergame == supergame[s], ]$round))
      
      for (p in seq_along(period)) {
        
        if (supergame[s] == 1 && period[p] == 1) {
          
          localdata[localdata$supergame == supergame[s] & localdata$round == period[p], "Output"] <- prior
          
        } else if (supergame[s] > 1 && period[p] == 1) {
          
          prev_supergame <- supergame[s - 1]
          max_round <- max(localdata[localdata$supergame == prev_supergame, ]$round)
          localdata[localdata$supergame == supergame[s] & localdata$round == period[p], "Output"] <- 
            ((lambda * localdata[localdata$supergame == prev_supergame & localdata$round == max_round, "Output"]) + ((1-lambda) * localdata[localdata$supergame == prev_supergame & localdata$round == max_round, varofinterest]))
          
        } else {
          
          localdata[localdata$supergame == supergame[s] & localdata$round == period[p], "Output"] <- 
            ((lambda * localdata[localdata$supergame == supergame[s] & localdata$round == period[p - 1], "Output"]) + ((1-lambda) * localdata[localdata$supergame == supergame[s] & localdata$round == period[p - 1], varofinterest]))
          
        }
        
        print(id[i])
        
      }
    }
    return(localdata)
  }
  return(result$Output)
}

#construct belief Markov Chain considering the variable of interest per supergame (PARALLEL)
belief_construction_markov_parallel_B2 <- function(metadata, lambda, varofinterest){
  
  metadata$Output <- NA
  prior <- 0.5
  id <- unique(metadata$id)
  
  # Parallel computation
  result <- foreach(i = 1:length(id), .combine = rbind, .packages = c("tidyverse")) %dopar% {
    
    localdata <- metadata[metadata$id == id[i], ]
    
    supergame <- sort(unique(localdata$supergame))
    
    for (s in seq_along(supergame)) {
      
      if (supergame[s] == 1) {
        
        localdata[localdata$supergame == supergame[s], "Output"] <- prior
        
      } else {
        
        supergame_of_interest <- supergame[s-1]
        prior_posterior <- unique(localdata[localdata$supergame == supergame_of_interest, "Output"])
        
        
        localdata[localdata$supergame == supergame[s], "Output"] <- 
          ((lambda * unique(localdata[localdata$supergame == supergame_of_interest, "Output"])) + ((1-lambda) * mean(localdata[localdata$supergame == supergame_of_interest, varofinterest])))
        
        
      }
      
    }
    return(localdata)
  }
  return(result$Output)
}

#construct belief Bayesian Updating considering the rounds (PARALLEL)
belief_construction_bayes_parallel_B1 <- function(metadata, beta, varofinterest){
  
  metadata$Output <- NA
  prior <- 0.5
  id <- unique(metadata$id)
  
  # Parallel computation
  result <- foreach(i = 1:length(id), .combine = rbind, .packages = c("tidyverse")) %dopar% {
    
    localdata <- metadata[metadata$id == id[i], ]
    
    supergame <- sort(unique(localdata$supergame))
    
    for (s in seq_along(supergame)) {
      
      period <- sort(unique(localdata[localdata$supergame == supergame[s], ]$round))
      
      for (p in seq_along(period)) {
        
        if (supergame[s] == 1 && period[p] == 1) {
          
          localdata[localdata$supergame == supergame[s] & localdata$round == period[p], "Output"] <- (1-prior)*(1-beta)+(prior*beta)
          
        } else {
          
          if(supergame[s] > 1 && period[p] == 1){
            
            supergame_of_interest <- supergame[s-1]
            round_of_interest <- max(localdata[localdata$supergame == supergame_of_interest, ]$round)
            
          }else{
            
            supergame_of_interest <- supergame[s]
            round_of_interest <-  period[p-1]
            
          }
          
          belief_t_minus_one <- localdata[localdata$supergame == supergame_of_interest & localdata$round == round_of_interest, "Output"]
          
          update_other_d <- ((1-belief_t_minus_one)*(1-beta))/((1-belief_t_minus_one)*(1-beta)+belief_t_minus_one*(beta))
          #update_other_d <- ((1-belief_t_minus_one)*beta)/((1-belief_t_minus_one)*beta+belief_t_minus_one*(1-beta))
          update_other_c <- (belief_t_minus_one*beta)/(belief_t_minus_one*beta+((1-belief_t_minus_one)*(1-beta)))
          
          observe_c <- update_other_d*(1-beta)+update_other_c*beta
          observe_d <- update_other_c*(1-beta)+update_other_d*beta
          
          c_t_minus_one <- localdata[localdata$supergame == supergame_of_interest & localdata$round == round_of_interest, varofinterest]
          
          localdata[localdata$supergame == supergame[s] & localdata$round == period[p], "Output"] <- observe_c*c_t_minus_one+observe_d*(1-c_t_minus_one)
          
          
          
        }
        
      }
    }
    return(localdata)
  }
  return(result$Output)
}

#construct belief gamma-Weighted Beliefs considering the rounds (PARALLEL)
belief_construction_weighted_parallel_B1 <- function(metadata, gamma, varofinterest){
  
  metadata$Output <- NA
  prior <- 0.5
  id <- unique(metadata$id)
  
  # Parallel computation
  result <- foreach(i = 1:length(id), .combine = rbind, .packages = c("tidyverse")) %dopar% {
    
    localdata <- metadata[metadata$id == id[i], ]
    
    supergame <- sort(unique(localdata$supergame))
    
    for (s in seq_along(supergame)) {
      
      period <- sort(unique(localdata[localdata$supergame == supergame[s], ]$round))
      
      for (p in seq_along(period)) {
        
        if (supergame[s] == 1 && period[p] == 1) {
          
          localdata[localdata$supergame == supergame[s] & localdata$round == period[p], "Output"] <- prior
          
        } else {
          
          t <- nrow(localdata[(localdata$supergame < supergame[s] | (localdata$supergame == supergame[s] & localdata$round < period[p])),])
          u <- 1:(t-1)
          
          #varofinterest in the previous round
          part_one <- localdata[(localdata$supergame < supergame[s] | (localdata$supergame == supergame[s] & localdata$round < period[p])),][t,varofinterest]
          
          part_two <- (sum((gamma)^(u)*(localdata[(localdata$supergame < supergame[s] | (localdata$supergame == supergame[s] & localdata$round < period[p])),][t-u,varofinterest])))
          
          nominator <- part_one + part_two
          
          denominator <- (1 + sum((gamma)^(u)))
          
          #prior in this round is
          localdata[localdata$supergame == supergame[s] & localdata$round == period[p], "Output"] <- (nominator/denominator)
          
        }
        
        print(id[i])
        
      }
    }
    return(localdata)
  }
  return(result$Output)
}

#construct belief gamma-Weighted Beliefs considering the variable of interest per supergame (PARALLEL)
belief_construction_weighted_parallel_B2 <- function(metadata, gamma, varofinterest){
  
  metadata$Output <- NA
  prior <- 0.5
  id <- unique(metadata$id)
  
  # Parallel computation
  result <- foreach(i = 1:length(id), .combine = rbind, .packages = c("tidyverse")) %dopar% {
    
    localdata <- metadata[metadata$id == id[i], ]
    
    supergame <- sort(unique(localdata$supergame))
    
    for (s in seq_along(supergame)) {

      if (supergame[s] == 1) {
        
        localdata[localdata$supergame == supergame[s], "Output"] <- prior
        
      } else {
        
        if (supergame[s] == 2) {
          
          localdata[localdata$supergame == supergame[s], "Output"] <- mean(localdata[localdata$supergame == 1,][,varofinterest])
          
        }else{

          t <- max(unique(localdata[localdata$supergame < supergame[s],][,"supergame"]))
          u <- 1:(t-1)
          
          #varofinterest in the previous supergame
          part_one <- mean(localdata[localdata$supergame == t,][,varofinterest])
          
          part_two <- (sum((gamma)^(u)*(mean(localdata[localdata$supergame == supergame[t-u],][,varofinterest]))))
          
          nominator <- part_one + part_two
          
          denominator <- (1 + sum((gamma)^(u)))
          
          #prior in this supergane is
          localdata[localdata$supergame == supergame[s], "Output"] <- (nominator/denominator)
          
        }
        
      }
      
      print(id[i])
        
    }
    return(localdata)
  }
  return(result$Output)
}

#construct belief Learning Model Beliefs considering the variable of interest per supergame (PARALLEL)
belief_construction_learning_parallel_B2 <- function(metadata, theta, varofinterest){
  
  metadata$Output <- NA
  metadata$Experience <- NA
  
  #the values for alpha and beta comes from Fudenberg and Gustav Karreskog Rehbinder (2024) in Predicting Cooperation using Learning Models
  alpha=0#-0.268
  beta=0#1.291
  #theta, i.e., lambda from the paper 0.182 

  id <- unique(metadata$id)
  
  # Parallel computation
  result <- foreach(i = 1:length(id), .combine = rbind, .packages = c("tidyverse")) %dopar% {
    
    localdata <- metadata[metadata$id == id[i], ]
    
    supergame <- sort(unique(localdata$supergame))
    
    for (s in seq_along(supergame)) {
      
      Delta_rd = localdata[localdata$supergame == supergame[s] & localdata$round == 1, "delta"] - localdata[localdata$supergame == supergame[s] & localdata$round == 1, "rd"]
      
      if (supergame[s] == 1) {
        
        e=0
        
        localdata[localdata$supergame == supergame[s], "Output"] <- (1)/(1+exp(-(alpha + beta*Delta_rd +e)))
        
        localdata[localdata$supergame == supergame[s], "Experience"] <- e
        
      } else {
        
        #a <- ifelse(localdata[localdata$supergame == supergame[s-1] & localdata$round == 1, varofinterest] == 1, 1, -1)
        a <- ifelse(mean(localdata[localdata$supergame == supergame[s-1], varofinterest]) >= 0.5, 1, -1)
        
        v <- sum(localdata[localdata$supergame == supergame[s-1], "Payoff"])
        
        e_before <- localdata[localdata$supergame == supergame[s-1] & localdata$round == 1, "Experience"]
        
        e <- (theta*a*v)+e_before
        
        localdata[localdata$supergame == supergame[s], "Experience"] <- e
        
        localdata[localdata$supergame == supergame[s], "Output"] <- (1)/(1+exp(-(alpha + beta*Delta_rd +e)))
        
      }
      
      print(id[i])
        

    }
    return(localdata)
  }
  return(result$Output)
}

#constructed beliefs
cdf_beliefs <- function(metadata, lambda, beliefofinterest, belieflevelofinterest, colour_learning){
  
  size=1
  
  cdf <- ggplot(data = metadata[metadata$supergame > 1 & metadata$round == 1,],) 
  
  for (i in 1:length(lambda)) {
    
    #plot cdf
    cdf <- cdf + stat_ecdf(aes_string(x=paste(beliefofinterest, lambda[i], belieflevelofinterest,sep="")),
                           geom = "step", size = size, color = colour_learning[i], pad = FALSE)

    
  }
  
  cdf<-cdf+
    theme_bw() +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    labs(y = "Cumulative Distribution Function", x="Belief") +
    theme(legend.position="bottom")
  
  cdf
  
  ggsave(paste("Figures Meta Data/Figure_ConstructedBeliefs_", beliefofinterest, belieflevelofinterest, ".pdf", sep = ""), width = 3.5, height = 3.5)
  
  
}

#create partial effect data
PartialEffectRegression <- function(metadata, lambda, ranges, beliefofinterest, belieflevelofinterest){
  
  PartialEffectData  <- data.frame()
  
  for (i in 1:length(lambda)){
    
    PartialEffectData[i,"LearningParameter"] <- lambda[i]
    
    for(j in 1:(length(ranges)-1)){
      
      #Model for g and l
      Model <- probitmfx(formula = Cooperation ~ g + l, 
                         data = metadata[metadata[paste(beliefofinterest, lambda[i], belieflevelofinterest,sep="")] > ranges[j] & metadata[paste(beliefofinterest, lambda[i], belieflevelofinterest,sep="")] <= ranges[j+1] & metadata$round == 1 & metadata$supergame > 1,],
                         atmean = FALSE, robust = TRUE, clustervar1 = "id")
      
      
      PartialEffectData[i,paste("M", j, "Effect","g",sep="")] <- Model$mfxest["g","dF/dx"]
      PartialEffectData[i,paste("M", j, "SE","g",sep="")] <- Model$mfxest["g","Std. Err."]
      PartialEffectData[i,paste("M", j, "z","g",sep="")] <- Model$mfxest["g","z"]
      PartialEffectData[i,paste("M", j, "p-value","g",sep="")] <- Model$mfxest["g","P>|z|"]
      
      PartialEffectData[i,paste("M", j, "Effect","l",sep="")] <- Model$mfxest["l","dF/dx"]
      PartialEffectData[i,paste("M", j, "SE","l",sep="")] <- Model$mfxest["l","Std. Err."]
      PartialEffectData[i,paste("M", j, "z","l",sep="")] <- Model$mfxest["l","z"]
      PartialEffectData[i,paste("M", j, "p-value","l",sep="")] <- Model$mfxest["l","P>|z|"]
      

      #Model for differences
      Model<- probitmfx(formula = Cooperation ~ g+gplusl, 
                         data = metadata[metadata[paste(beliefofinterest, lambda[i], belieflevelofinterest,sep="")] > ranges[j] & metadata[paste(beliefofinterest, lambda[i], belieflevelofinterest,sep="")] <= ranges[j+1] & metadata$round == 1 & metadata$supergame > 1,],
                         atmean = FALSE, robust = TRUE, clustervar1 = "id")
      
      PartialEffectData[i,paste("M", j, "Effect","diff",sep="")] <- Model$mfxest["g","dF/dx"]
      PartialEffectData[i,paste("M", j, "SE","diff",sep="")] <- Model$mfxest["g","Std. Err."]
      PartialEffectData[i,paste("M", j, "z","diff",sep="")] <- Model$mfxest["g","z"]
      PartialEffectData[i,paste("M", j, "p-value","diff",sep="")] <- Model$mfxest["g","P>|z|"]
      
      
    }
    
  }
  
  # PartialEffectData[length(lambda)+1,"LearningParameter"] <- 1
  # 
  # for (j in 1:(length(ranges)-1)) {
  #   
  #   PartialEffectData[length(lambda)+1, paste("M", j, "Effect","g",sep="")] <- mean(PartialEffectData[1:length(lambda),paste("M", j, "Effect","g",sep="")])
  #   PartialEffectData[length(lambda)+1, paste("M", j, "SE","g",sep="")] <- mean(PartialEffectData[1:length(lambda),paste("M", j, "SE","g",sep="")])
  #   
  #   PartialEffectData[length(lambda)+1, paste("M", j, "Effect","l",sep="")] <- mean(PartialEffectData[1:length(lambda),paste("M", j, "Effect","l",sep="")])
  #   PartialEffectData[length(lambda)+1, paste("M", j, "SE","l",sep="")] <- mean(PartialEffectData[1:length(lambda),paste("M", j, "SE","l",sep="")])
  #   
  #   PartialEffectData[length(lambda)+1, paste("M", j, "Effect","diff",sep="")] <- mean(PartialEffectData[1:length(lambda),paste("M", j, "Effect","diff",sep="")])
  #   PartialEffectData[length(lambda)+1, paste("M", j, "SE","diff",sep="")] <- mean(PartialEffectData[1:length(lambda),paste("M", j, "SE","diff",sep="")])
  #   
  # }
  
  return(PartialEffectData)
  
}

#create partial effect plot of gain layer
plotLayers<-function(PartialEffectData, ranges, layer, colour_learning, y_axis_title){
  
  breaks <- c(0.1,0.3,0.5,0.7,0.9)
  
  p<-ggplot()+geom_hline(yintercept = 0, linetype = "dashed")
  
  for(i in 1:length(PartialEffectData[,"LearningParameter"])){
    
    for(j in 1:(length(ranges)-1)){
      
      estimate <- PartialEffectData[i, paste("M", j, "Effect", layer, sep="")]
      
      se <- PartialEffectData[i, paste("M", j, "SE", layer,sep="")]
      
      nudge <- 0.05
      
      if(i == length(PartialEffectData[,"LearningParameter"])){
        
        shift <- 0
        
      }else{
        
        shift <- ((i-mean(c(1:(length(PartialEffectData[,"LearningParameter"])-1))))*nudge)
        
      }
      
      position <- breaks[j] + shift
      
      
      shape <- 19
      
      size_point = 3
      
      size_line = 1
      
      p<-p+
        annotate("point", x = position, y = estimate, color = colour_learning[i], fill = colour_learning[i],shape=shape, size=size_point) +
        annotate("errorbar", 
                 x = position, 
                 y = estimate, 
                 ymin = estimate - se*1.96, 
                 ymax = estimate + se*1.96 , 
                 color = colour_learning[i], size = size_line, width = 0.03)
      
      
      if(i == length(PartialEffectData[,"LearningParameter"]) & j > 1){
        
        p<-p+
          annotate("segment", x = breaks[j-1]+shift, xend = breaks[j]+shift,
                   y = PartialEffectData[i, paste("M", j-1, "Effect", layer,sep="")], yend = PartialEffectData[i, paste("M", j, "Effect", layer,sep="")],
                   color = colour_learning[i], size = size_line)
        
      }else{
        
      }
      
      
    }
    
    
  }
  
  p<-p+
    theme_bw() +
    coord_cartesian(xlim = c(0,1), ylim = c(-0.3,0.3)) +
    scale_x_continuous(breaks = breaks) +#scale_x_continuous(name = NULL, breaks = breaks) +
    xlab("Belief") + 
    ylab(y_axis_title)+
    theme(legend.position="bottom")
  
  return(p)
  
}

#calculate partial effects
PartialEffect <- function(metadata, lambda, ranges, beliefofinterest, belieflevelofinterest){
  
  #create data
  PartialEffectData <- PartialEffectRegression(metadata, lambda, ranges, beliefofinterest, belieflevelofinterest)
  
  #plot layers
  plotLayers(PartialEffectData, ranges, layer = "g", colour_learning = c("#B3B3B3","#B3B3B3","#d95f02"), y_axis_title ="Partial Effect of the Gain on Cooperation")
  ggsave(paste("Figures Meta Data/Figure_PartialEffects_", "g","_", beliefofinterest, belieflevelofinterest, ".pdf", sep = ""), width = 2.5, height = 3.5)
  
  plotLayers(PartialEffectData, ranges, layer = "l", colour_learning = c("#B3B3B3","#B3B3B3","#7570b3"), y_axis_title ="Partial Effect of the Loss on Cooperation")
  ggsave(paste("Figures Meta Data/Figure_PartialEffects_", "l","_", beliefofinterest, belieflevelofinterest, ".pdf", sep = ""), width = 2.5, height = 3.5)
  
  plotLayers(PartialEffectData, ranges, layer = "diff", colour_learning = c("#B3B3B3","#B3B3B3","#1b9e77"), y_axis_title ="Partial Effect of the Gain - Partial Effect of the Loss")
  ggsave(paste("Figures Meta Data/Figure_PartialEffects_", "diff","_", beliefofinterest, belieflevelofinterest, ".pdf", sep = ""), width = 2.5, height = 3.5)
  

  return(PartialEffectData)
  
}

#average cooperation per constructed belief
CooperationBeliefAverage <- function(metadata, lambda, ranges, beliefofinterest, belieflevelofinterest){
  
  CooperationBeliefData  <- data.frame()
  
  for (i in 1:length(lambda)){
    
    CooperationBeliefData[i,"LearningParameter"] <- lambda[i]
    
    for(j in 1:(length(ranges)-1)){
      
      #Model for mean and se belief
      Model <- metadata[metadata[paste(beliefofinterest, lambda[i], belieflevelofinterest,sep="")] > ranges[j] & metadata[paste(beliefofinterest, lambda[i], belieflevelofinterest,sep="")] <= ranges[j+1] & metadata$round == 1 & metadata$supergame > 1,]%>%
        group_by(id)%>%
        summarise(mean = mean(Markov0.25B1), sd = sd(Markov0.25B1))
      
      CooperationBeliefData[i,paste("M", j, "Mean",sep="")] <- mean(Model$mean)
      CooperationBeliefData[i,paste("M", j, "SE",sep="")] <- (sd(Model$mean)/sqrt(length(Model$mean)))
      
      
    }
    
  }
  
  return(CooperationBeliefData)
}

#plot Cooperation Belief Average
plotCooperationBelief <- function(CooperationBeliefData, ranges, colour_learning){
  
  breaks <- c(0.1,0.3,0.5,0.7,0.9)
  
  
  p<-ggplot()
  
  for(i in 1:length(CooperationBeliefData[,"LearningParameter"])){
    
    for(j in 1:(length(ranges)-1)){
      
      estimate <- CooperationBeliefData[i, paste("M", j, "Mean", sep="")]
      
      se <- CooperationBeliefData[i, paste("M", j, "SE",sep="")]
      
      nudge <- 0.05
      
      if(i == length(CooperationBeliefData[,"LearningParameter"])){
        
        shift <- 0
        
      }else{
        
        shift <- ((i-mean(c(1:(length(CooperationBeliefData[,"LearningParameter"])-1))))*nudge)
        
      }
      
      position <- breaks[j] + shift
      
      
      shape <- 19
      
      size_point = 3
      
      size_line = 1
      
      p<-p+
        annotate("point", x = position, y = estimate, color = colour_learning[i], fill = colour_learning[i],shape=shape, size=size_point) +
        annotate("errorbar", 
                 x = position, 
                 y = estimate, 
                 ymin = estimate - se*1.96, 
                 ymax = estimate + se*1.96 , 
                 color = colour_learning[i], size = size_line, width = 0.03)
      
      
      if(i == length(CooperationBeliefData[,"LearningParameter"]) & j > 1){
        
        p<-p+
          annotate("segment", x = breaks[j-1]+shift, xend = breaks[j]+shift,
                   y = CooperationBeliefData[i, paste("M", j-1, "Mean",sep="")], yend = CooperationBeliefData[i, paste("M", j, "Mean",sep="")],
                   color = colour_learning[i], size = size_line)
        
      }else{
        
      }
      
      
    }
    
    
  }
  
  p<-p+
    theme_bw() +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    scale_x_continuous(name = NULL, breaks = breaks) +
    xlab("Belief") + 
    ylab("Average Cooperation")+
    theme(legend.position="bottom")
  
  
  return(p)
  
}

#calculate cooperation per constructed belief
CooperationBelief <- function(metadata, lambda, ranges, beliefofinterest, belieflevelofinterest){
  
  #create cooperation belief data
  CooperationBeliefData <- CooperationBeliefAverage(metadata, lambda, ranges, beliefofinterest, belieflevelofinterest)
  
  #plot the cooperation belief data
  plotCooperationBelief(CooperationBeliefData, ranges, colour_learning = c("#B3B3B3","#B3B3B3","#377eb8"))
  ggsave(paste("Figures Meta Data/Figure_CooperationBeliefs_", beliefofinterest, belieflevelofinterest, ".pdf", sep = ""), width = 3.5, height = 3.5)
  
  
}



####################
#
# Main
#
#######################

#load data
metadata <- load_data(wd, filename)


# belief construction
ranges <- c(-1,0.2,0.4,0.6,0.8,1)

lambda <- c(0.25, 0.75, 0.5)
beta <- c(0.95,0.99,0.975)
gamma <- c(0.25,0.75, 0.5)
theta <- c(0.03,0.07,0.05)



cl = makeCluster(14)
registerDoParallel(cl)
for (i in 1:length(lambda)){
  metadata[paste("Markov", lambda[i], "B1",sep="")] <- belief_construction_markov_parallel_B1(metadata, lambda = lambda[i], varofinterest = "CooperationRate")
}
for (i in 1:length(beta)){
  metadata[paste("Bayes", beta[i], "B1",sep="")] <- belief_construction_bayes_parallel_B1(metadata, beta = beta[i], varofinterest = "MutualCooperation")
}
for (i in 1:length(gamma)){
  metadata[paste("Weighted", gamma[i], "B1",sep="")] <- belief_construction_weighted_parallel_B1(metadata, gamma = gamma[i], varofinterest = "CooperationRate")
}
for (i in 1:length(theta)){
  metadata[paste("Learning", theta[i], "B2",sep="")] <- belief_construction_learning_parallel_B2(metadata, theta = theta[i], varofinterest = "Cooperation")
}
stopCluster(cl)

summary(lm(Cooperation ~ Markov0.5B1 + g + Markov0.5B1*g, data = metadata[metadata$supergame > 1 & metadata$round == 1,]))


#constructed beliefs
cdf_beliefs(metadata = metadata, lambda = lambda, beliefofinterest = "Markov", belieflevelofinterest = "B1", colour_learning = c("#B3B3B3","#B3B3B3","#000000"))
cdf_beliefs(metadata = metadata, lambda = beta, beliefofinterest = "Bayes", belieflevelofinterest = "B1", colour_learning = c("#B3B3B3","#B3B3B3","#000000"))
cdf_beliefs(metadata = metadata, lambda = gamma, beliefofinterest = "Weighted", belieflevelofinterest = "B1", colour_learning = c("#B3B3B3","#B3B3B3","#000000"))
cdf_beliefs(metadata = metadata, lambda = theta, beliefofinterest = "Learning", belieflevelofinterest = "B2", colour_learning = c("#B3B3B3","#B3B3B3","#000000"))

#average cooperation per belief
CooperationBelief(metadata = metadata, lambda = lambda, ranges = ranges, beliefofinterest = "Markov", belieflevelofinterest = "B1")
CooperationBelief(metadata = metadata, lambda = beta, ranges = ranges, beliefofinterest = "Bayes", belieflevelofinterest = "B1")
CooperationBelief(metadata = metadata, lambda = gamma, ranges = ranges, beliefofinterest = "Weighted", belieflevelofinterest = "B1")
CooperationBelief(metadata = metadata, lambda = theta, ranges = ranges, beliefofinterest = "Learning", belieflevelofinterest = "B2")

#partial effect 
PartialEffect(metadata = metadata, lambda = lambda, ranges = ranges, beliefofinterest = "Markov", belieflevelofinterest = "B1")
PartialEffect(metadata = metadata, lambda = beta, ranges = ranges, beliefofinterest = "Bayes", belieflevelofinterest = "B1")
PartialEffect(metadata = metadata, lambda = gamma, ranges = ranges, beliefofinterest = "Weighted", belieflevelofinterest = "B1")
PartialEffect(metadata = metadata, lambda = theta, ranges = ranges, beliefofinterest = "Learning", belieflevelofinterest = "B2")





#






