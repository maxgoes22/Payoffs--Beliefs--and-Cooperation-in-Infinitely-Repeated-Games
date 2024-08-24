
##################################
#
# This script covers the function for the Meta Study
#
##################################

#load the data
load_data <- function(wd, filename){
  
  #read meta data
  metadata <- read.delim(paste(wd, filename, sep=""), header = TRUE, sep = "\t", dec = ".")
  
  #remove one-shot games
  metadata <- metadata[metadata$delta > 0,]

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

#plot the learning path
plot_learning_path <- function(lambda){
  
  #install.packages("ggplot2")
  library("ggplot2")
  
  nperiod = 50
  sizebad = 0.5
  startingbelief = 0.5
  errorrate = 0.1
  lowtype = 0.2
  hightype = 0.8
  nsim = 100
  
  
  outerresult = matrix(, nrow = nperiod*2, ncol = nsim)
  
  for(k in seq(1, nsim, 1)){
    
    set.seed(k+1); errori = ifelse(runif(nperiod*2, 0, 1) >= (1-errorrate), 1, 0) 
    set.seed(k+2); coopj = runif(nperiod*2, 0, 1) 
    
    beliefdata = data.frame(
      t = rep(seq(1,nperiod,1), times = 2),
      othertype = rep(c(lowtype,hightype), each = nperiod),
      pi = -999,
      ai = -999,
      aj = -999, 
      c = -999,
      lambda = lambda,
      errori = errori,
      coopj = coopj
    )
    
    for(i in 1:nrow(beliefdata)){
      
      if(beliefdata[i,"othertype"] == lowtype){
        t = i
      }else if (beliefdata[i,"othertype"] == hightype){
        t = i - nperiod
      }
      
      # Player i's belief and action
      beliefdata[i,"pi"] = if(t == 1){startingbelief}else{(beliefdata[i,"lambda"]*beliefdata[i-1,"pi"])+((1-beliefdata[i,"lambda"])*beliefdata[i-1,"c"])}
      beliefdata[i,"ai"] = if(beliefdata[i,"pi"] >= sizebad){1}else{0}
      beliefdata[i,"ai"] = if(beliefdata[i,"errori"] == 1){1-beliefdata[i,"ai"]}else{beliefdata[i,"ai"]}
      
      # Player j's action using cooperation probability
      beliefdata[i,"aj"] = if(beliefdata[i,"othertype"] > beliefdata[i,"coopj"]){1}else{0}
      
      # Cooperation rate
      beliefdata[i,"c"] = (beliefdata[i,"ai"] + beliefdata[i,"aj"])/2
    }
    
    outerresult[, k] = beliefdata$pi
    
  }
  
  
  # Plot
  rowVar <- function(x, ...) {
    (rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1))
  }
  
  beliefresults = data.frame(
    t = rep(seq(1,nperiod,1), times = 2),
    othertype = as.factor(rep(c(lowtype,hightype), each = nperiod)),
    lambda = lambda,
    pi = rowMeans(outerresult),
    piVar = rowVar(outerresult)
  )
  
  p <- ggplot() +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_ribbon(data = beliefresults, aes(x = t, ymin = pi-piVar, ymax = pi + piVar, fill = othertype), alpha = 0.5) +
    geom_line(data = beliefresults, aes(x = t, y = pi, color = othertype), linewidth = 1) +
    theme_bw() +
    coord_cartesian(xlim = c(0, nperiod), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, nperiod, nperiod / 5), name = "Period") +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Player i's Belief") +
    scale_color_discrete(name = "j's P(C)") +
    scale_fill_discrete(guide = "none") +
    theme(legend.position = "bottom")
  
  return(p)
  
}

#construct belief Markov Chain considering the rounds (PARALLEL)
belief_construction_markov_parallel <- function(metadata, lambda, varofinterest){
  
  #install.packages("foreach")
  library("foreach")
  #install.packages("tidyverse")
  library("tidyverse")
  
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

#create histogram plot of the constructed beliefs
plot_histogram <- function(metadata, lambda, appendix = FALSE){
  
  #install.packages("rlang")
  library(rlang)
  #install.packages("latex2exp")
  library("latex2exp")
  #install.packages("ggplot2")
  library("ggplot2")
  
  x_variable <- paste("Markov",lambda,sep="")
  
  x_title <- paste("Markov Belief with $\\lambda = ",lambda,"$",sep="")
  
  breaks_y <- if(appendix == TRUE){seq(0,10000,2000)}else{seq(0,15000,2000)}
  
  y_lim <- if(appendix == TRUE){c(0,12000)}else{c(0,15000)}
  
  # Histogram
  p <- ggplot(data = metadata[metadata$round == 1, ], aes(x = !! sym(x_variable))) +
    geom_histogram(color = "#377eb8", fill = "#377eb8", linewidth = 0.5, alpha = 0.25, bins = 10) +
    theme_bw() +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim=y_lim) +
    scale_x_continuous(name = TeX(x_title), breaks = c(0,0.25,0.5,0.75,1)) +
    scale_y_continuous(breaks = breaks_y, name = "Count")
  
  
  return(p)
  
}

#create confusion matrix plot of the constructed beliefs and print the accuracy
plot_confusion_matrix <- function(metadata, lambda, appendix = FALSE){
  
  #install.packages("latex2exp")
  library("latex2exp")
  #install.packages("ggplot2")
  library("ggplot2")
  
  var_of_interest <- paste("Markov",lambda,sep='')
  
  #calculate the critical discount factor 
  metadata$MarkovCriticalDelta <- (metadata[,var_of_interest] * (metadata$g - metadata$l) + metadata$l) / (metadata[,var_of_interest] * (1 + metadata$g - metadata$l) + metadata$l)
  
  #calculate equilibrium selection
  metadata$MarkovCriticalDeltaCoop <- ifelse(metadata$delta >= metadata$MarkovCriticalDelta, 1, 0)
  
  #calculate confusion
  MarkovConfusionFreq <- c(
    nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 0 & metadata$Cooperation == 0, ]) / nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 0, ]),
    nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 1 & metadata$Cooperation == 0, ]) / nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 1, ]),
    nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 0 & metadata$Cooperation == 1, ]) / nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 0, ]),
    nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 1 & metadata$Cooperation == 1, ]) / nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == 1, ])
  )
  
  #create confusion data
  MarkovConfusionData <- data.frame(
    PredCoop = factor(c(0,1,0,1)),
    Cooperation = factor(c(0,0,1,1)),
    Freq = round(MarkovConfusionFreq,2)
  )
  
  #print the accuracy
  accuracy <- round(nrow(metadata[metadata$round == 1 & metadata$MarkovCriticalDeltaCoop == metadata$Cooperation, ]) / nrow(metadata[metadata$round == 1, ]), 2)
  print(paste("The accuracy of the constructed Markov belief with lambda of ",lambda," equals: ",accuracy,sep=""))
  
  #plot the confusion matrix
  legend <- ifelse(appendix == TRUE, "off", "bottom")
  x_title <- paste("$\\delta \\geq \\delta^*(p = \\mu_{\\lambda = ",lambda,"},g,l)$", sep='')
  
  p <- ggplot(data = MarkovConfusionData, aes(x = PredCoop, y = Cooperation, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#377eb8", limits = c(0, 1),
                        name = "Col. Freq.", breaks = c(0, 0.5, 1),) +
    geom_text(aes(label = round(Freq, 3)), color = "black", size = 3) +
    scale_x_discrete(name = TeX(x_title)) +
    theme_bw() +
    theme(legend.position = legend, ) +
    coord_fixed() 
  
  
  return(p)
  
  
}

#create partial effect data
partial_effect_regression <- function(metadata, lambda, ranges, beliefofinterest){
  
  #install.packages("mfx")
  library("mfx")
  
  PartialEffectData  <- data.frame()
  
  for (i in 1:length(lambda)){
    
    PartialEffectData[i,"LearningParameter"] <- lambda[i]
    
    for(j in 1:(length(ranges)-1)){
      
      #Model for g and l
      Model <- probitmfx(formula = Cooperation ~ g + l, 
                         data = metadata[metadata[paste(beliefofinterest, lambda[i],sep="")] > ranges[j] & metadata[paste(beliefofinterest, lambda[i],sep="")] <= ranges[j+1],],
                         atmean = FALSE, robust = TRUE, clustervar1 = "id")
      
      
      PartialEffectData[i,paste("M", j, "Effect","g",sep="")] <- Model$mfxest["g","dF/dx"]
      PartialEffectData[i,paste("M", j, "SE","g",sep="")] <- Model$mfxest["g","Std. Err."]
      PartialEffectData[i,paste("M", j, "z","g",sep="")] <- Model$mfxest["g","z"]
      PartialEffectData[i,paste("M", j, "p-value","g",sep="")] <- Model$mfxest["g","P>|z|"]
      PartialEffectData[i,paste("M", j, "N","g",sep="")] <- nobs(Model$fit)
      
      PartialEffectData[i,paste("M", j, "Effect","l",sep="")] <- Model$mfxest["l","dF/dx"]
      PartialEffectData[i,paste("M", j, "SE","l",sep="")] <- Model$mfxest["l","Std. Err."]
      PartialEffectData[i,paste("M", j, "z","l",sep="")] <- Model$mfxest["l","z"]
      PartialEffectData[i,paste("M", j, "p-value","l",sep="")] <- Model$mfxest["l","P>|z|"]
      PartialEffectData[i,paste("M", j, "N","l",sep="")] <- nobs(Model$fit)
      

      #Model for differences
      Model<- probitmfx(formula = Cooperation ~ g+gplusl, 
                         data = metadata[metadata[paste(beliefofinterest, lambda[i],sep="")] > ranges[j] & metadata[paste(beliefofinterest, lambda[i],sep="")] <= ranges[j+1],],
                         atmean = FALSE, robust = TRUE, clustervar1 = "id")
      
      PartialEffectData[i,paste("M", j, "Effect","diff",sep="")] <- Model$mfxest["g","dF/dx"]
      PartialEffectData[i,paste("M", j, "SE","diff",sep="")] <- Model$mfxest["g","Std. Err."]
      PartialEffectData[i,paste("M", j, "z","diff",sep="")] <- Model$mfxest["g","z"]
      PartialEffectData[i,paste("M", j, "p-value","diff",sep="")] <- Model$mfxest["g","P>|z|"]
      PartialEffectData[i,paste("M", j, "N","diff",sep="")] <- nobs(Model$fit)
      
      PartialEffectData[i,paste("M", j, "Effect","gplusl",sep="")] <- Model$mfxest["gplusl","dF/dx"]
      PartialEffectData[i,paste("M", j, "SE","gplusl",sep="")] <- Model$mfxest["gplusl","Std. Err."]
      PartialEffectData[i,paste("M", j, "z","gplusl",sep="")] <- Model$mfxest["gplusl","z"]
      PartialEffectData[i,paste("M", j, "p-value","gplusl",sep="")] <- Model$mfxest["gplusl","P>|z|"]
      PartialEffectData[i,paste("M", j, "N","gplusl",sep="")] <- nobs(Model$fit)
      
    }
    
  }
  
  return(PartialEffectData)
  
}


# output partial effects regression data as table
table_regression_data <- function(regressiondata, lambda, asLaTeX){
  
  PartialEffectsTable = data.frame(
    Variable = c("Effectg", "SEg", "p-valueg", "Effectl", "SEl", "p-valuel", "Nl" ),
    Col1 = NA,
    Col2 = NA,
    Col3 = NA,
    Col4 = NA,
    Col5 = NA
  )
  
  for(var in c("Effectg", "SEg", "p-valueg", "Effectl", "SEl", "p-valuel", "Nl" )){
    for(M in 1:5){
      value <- regressiondata[regressiondata$LearningParameter == lambda, paste("M", M, var, sep = "")]
      PartialEffectsTable[PartialEffectsTable$Variable == var, M + 1] <- round(value, 4)
    }
  }
  
  PartialEffectsDiffTable = data.frame(
    Variable = c("Effectdiff", "SEdiff", "p-valuediff", "Effectgplusl", "SEgplusl", "p-valuegplusl", "Ngplusl" ),
    Col1 = NA,
    Col2 = NA,
    Col3 = NA,
    Col4 = NA,
    Col5 = NA
  )
  
  for(var in c("Effectdiff", "SEdiff", "p-valuediff", "Effectgplusl", "SEgplusl", "p-valuegplusl", "Ngplusl" )){
    for(M in 1:5){
      value <- regressiondata[regressiondata$LearningParameter == lambda, paste("M", M, var, sep = "")]
      PartialEffectsDiffTable[PartialEffectsDiffTable$Variable == var, M + 1] <- round(value, 4)
    }
  }
  
  if(asLaTeX){
    #install.packages("xtable")
    library("xtable")
    print("------------------M1------------------")
    print("Table of the average partial effects of gain and loss")
    print(xtable(PartialEffectsTable, digits = 4))
    print("------------------M2------------------")
    print("Table of the difference in the average partial effects of gain and loss")
    print(xtable(PartialEffectsDiffTable, digits = 4))
  }else{
    print("------------------M1------------------")
    print("Table of the average partial effects of gain and loss")
    print(PartialEffectsTable)
    print("------------------M2------------------")
    print("Table of the difference in the average partial effects of gain and loss")
    print(PartialEffectsDiffTable)
  }
  
}


#create partial effect plot of gain layer
plot_layers<-function(PartialEffectData, ranges, layer, colour_learning, y_axis_title, yrange){
  
  #install.packages("ggplot2")
  library("ggplot2")
  
  breaks <- c(0.1,0.3,0.5,0.7,0.9)
  
  p<-ggplot()+geom_hline(yintercept = 0, linetype = "dashed")
  
  for(i in 1:length(PartialEffectData[,"LearningParameter"])){
    
    for(j in 1:(length(ranges)-1)){
      
      estimate <- PartialEffectData[i, paste("M", j, "Effect", layer, sep="")]
      
      se <- PartialEffectData[i, paste("M", j, "SE", layer,sep="")]
      
      nudge <- 0.025
      
      if(i == length(PartialEffectData[,"LearningParameter"])){
        
        shift <- 0
        
      }else{
        
        if(i < mean(c(1:(length(PartialEffectData[,"LearningParameter"]))))){
          shift <- ((i-mean(c(1:(length(PartialEffectData[,"LearningParameter"]))))) * nudge)
        }else{
          shift <- ((i-mean(c(1:(length(PartialEffectData[,"LearningParameter"]))))) * nudge) + nudge
        }
        
      }
      
      position <- breaks[j] + shift
      
      
      shape <- 19
      
      size_point = if(length(PartialEffectData[,"LearningParameter"]) == 1){3}else{1}
      
      size_line = if(length(PartialEffectData[,"LearningParameter"]) == 1){1}else{0.5}
      
      p<-p+
        annotate("point", x = position, y = estimate, color = colour_learning[i], fill = colour_learning[i],shape=shape, size=size_point) +
        annotate("errorbar", 
                 x = position, 
                 y = estimate, 
                 ymin = estimate - se*1.96, 
                 ymax = estimate + se*1.96 , 
                 color = colour_learning[i], linewidth = size_line, width = 0.03)
      
      
      if(i == length(PartialEffectData[,"LearningParameter"]) & j > 1){
        
        p<-p+
          annotate("segment", x = breaks[j-1]+shift, xend = breaks[j]+shift,
                   y = PartialEffectData[i, paste("M", j-1, "Effect", layer,sep="")], yend = PartialEffectData[i, paste("M", j, "Effect", layer,sep="")],
                   color = colour_learning[i], linewidth = size_line)
        
      }else{
        
      }
      
      
    }
    
    
  }
  
  p<-p+
    theme_bw() +
    coord_cartesian(xlim = c(0,1), ylim = c(-yrange,yrange)) +
    scale_x_continuous(breaks = breaks) +
    xlab("Belief") + 
    ylab(y_axis_title)+
    theme(legend.position="bottom")
  
  return(p)
  
}


