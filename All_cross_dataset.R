rm(list = ls())
library(class)
library(FNN)
library(gtools)
library(R.matlab)
library(tseries)
library(fmsb)
# Adding or removing intercept does not effect eigen analysis

# Analysis <- readline(cat("What type of features? ('HOM','Eigen')\n"))
Analysis <- 'HOM'
DataSet_train <- readline(cat("What data set to train on? ('nback', 'gradCPT', 'ANT', 'HCP')\n"))
DataSet_test <- readline(cat("What data set to test on? ('nback', 'gradCPT', 'ANT', 'HCP')\n"))

options(digits=2)
dir.create("/Users/Mehraveh/Desktop/Final/",showWarnings = FALSE)
dir.create("/Users/Mehraveh/Desktop/Final/cross_datasets",showWarnings = FALSE)

if (DataSet_train == 'nback' || DataSet_test == 'nback')
{
  nback <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/nback_HOMs_5_ZTransformed.csv", header=FALSE)
  # nback <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/nback_HOMs_WeightedFC_5_ZTransformed.csv", header=FALSE)
  str(nback)
  colnames(nback)[1] <- "individual"
  colnames(nback)[12] <- c("dprime") #, "dp_2bk", "dp_Deg", "dp_avg")
  
  nback$TC <- rowMeans(cbind(nback$V2,nback$V3,nback$V6))
  nback$RTC <- rowMeans(cbind(nback$V4,nback$V5,nback$V7,nback$V8,nback$V9,nback$V10))
  nback$RC <- rowMeans(cbind(nback$V11))
  
  runs <- 5
  l <- nrow(nback)
  HOM <- readMat("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/EigenHOM_nback_5.mat", header=FALSE)
  HOM <- HOM$HOM
  
  Eignes <-matrix(0,l,runs)
  for(subj in 1:l)
  {
    tmp <- eigen(HOM[subj,,])
    Eignes[subj,] <- tmp$values
  }
  colnames(Eignes) <- c("eig1","eig2","eig3","eig4","eig5")
  nback <- cbind(nback, Eignes)
  x.nback <- nback[,c(2:11,13:ncol(nback))]
  y.nback <- nback$dprime
  nback_original<-nback
}


if (DataSet_train == 'gradCPT' || DataSet_test == 'gradCPT')
{
  gradCPT <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/gradCPT_HOMs_5_ZTransformed.csv", header=FALSE)
  # gradCPT <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/gradCPT_HOMs_WeightedFC_5_ZTransformed.csv", header=FALSE)
  str(gradCPT)
  colnames(gradCPT)[1] <- "individual"
  colnames(gradCPT)[12] <- "dprime"
  
  gradCPT$TC <- rowMeans(cbind(gradCPT$V2,gradCPT$V3,gradCPT$V6))
  gradCPT$RTC <- rowMeans(cbind(gradCPT$V4,gradCPT$V5,gradCPT$V7,gradCPT$V8,gradCPT$V9,gradCPT$V10))
  gradCPT$RC <- rowMeans(cbind(gradCPT$V11))
  
  runs <- 5
  l <- nrow(gradCPT)
  HOM <- readMat("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/EigenHOM_gradCPT_5.mat", header=FALSE)
  HOM <- HOM$HOM
  
  Eignes <-matrix(0,l,runs)
  for(subj in 1:l)
  {
    tmp <- eigen(HOM[subj,,])
    Eignes[subj,] <- tmp$values
  }
  colnames(Eignes) <- c("eig1","eig2","eig3","eig4","eig5")
  gradCPT <- cbind(gradCPT, Eignes)
  
  
  x.gradCPT <- gradCPT[,c(2:11,13:ncol(gradCPT))]
  y.gradCPT <- gradCPT$dprime
  gradCPT_original<-gradCPT
}


if (DataSet_train == 'ANT' || DataSet_test == 'ANT')
{
  runs <- 7
  if (runs == 5)
  {
    ANT <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/ANT_HOMs_5_ZTransformed.csv", header=FALSE)
    # ANT <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/ANT_HOMs_WeightedFC_5_ZTransformed.csv", header=FALSE)
    str(ANT)
    colnames(ANT)[1] <- "individual"
    
    # 5x5 matrix
    colnames(ANT)[12] <- "Mean_RT"
    colnames(ANT)[13] <- "Std_RT"
    colnames(ANT)[14] <- "Alert"
    colnames(ANT)[15] <- "Orient"
    colnames(ANT)[16] <- "Conflict"
    ANT$TC <- rowMeans(cbind(ANT$V2,ANT$V3,ANT$V6))
    ANT$RTC <- rowMeans(cbind(ANT$V4,ANT$V5,ANT$V7,ANT$V8,ANT$V9,ANT$V10))
    ANT$RC <- rowMeans(cbind(ANT$V11))
    
    l <- nrow(ANT)
    HOM <- readMat("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/EigenHOM_ANT_5.mat", header=FALSE)
    HOM <- HOM$HOM
    
    Eignes <-matrix(0,l,runs)
    for(subj in 1:l)
    {
      tmp <- eigen(HOM[subj,,])
      Eignes[subj,] <- tmp$values
    }
    colnames(Eignes) <- c("eig1","eig2","eig3","eig4","eig5")
    ANT <- cbind(ANT, Eignes)
    x.ANT <- ANT[,c(2:11,17:ncol(ANT))]
  }
  
  if (runs == 7)
  {
    ANT <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/ANT_HOMs_7_ZTransformed.csv", header=FALSE)
    # ANT <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/ANT_HOMs_WeightedFC_7_ZTransformed.csv", header=FALSE)
    str(ANT)
    colnames(ANT)[1] <- "individual"
    
    ### 7x7 matrix
    colnames(ANT)[23] <- "Mean_RT"
    colnames(ANT)[24] <- "Std_RT"
    colnames(ANT)[25] <- "Alert"
    colnames(ANT)[26] <- "Orient"
    colnames(ANT)[27] <- "Conflict"
    ANT$TC <- rowMeans(cbind(ANT$V2,ANT$V3,ANT$V4,ANT$V5,ANT$V8,ANT$V9,ANT$V10,ANT$V13,ANT$V14,ANT$V17))
    ANT$RTC <- rowMeans(cbind(ANT$V6,ANT$V7,ANT$V11,ANT$V12,ANT$V15,ANT$V16,ANT$V18,ANT$V19,ANT$V20,ANT$V21))
    ANT$RC <- rowMeans(cbind(ANT$V22))
    
    l <- nrow(ANT)
    HOM <- readMat("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/EigenHOM_ANT_7.mat", header=FALSE)
    HOM <- HOM$HOM
    
    Eignes <-matrix(0,l,runs)
    for(subj in 1:l)
    {
      tmp <- eigen(HOM[subj,,])
      Eignes[subj,] <- tmp$values
    }
    colnames(Eignes) <- c("eig1","eig2","eig3","eig4","eig5","eig6","eig7")
    ANT <- cbind(ANT, Eignes)
    x.ANT <- ANT[,c(2:22,28:ncol(ANT))]
  }
  
  
  # ANTbehav <- readline(cat("What ANT behavior? (Standard deviation of RT: 'std', RT coefficient of variation 'cv')\n"))
  ANTbehav <- 'cv'
  if (ANTbehav == 'cv')
  {
    ANT$Std_RT<-ANT$Std_RT/ANT$Mean_RT
  }
  
  if (ANTbehav == 'std')
  {
    behav_label = "Std RT"
  } else if (ANTbehav == 'cv') {
    behav_label = "RT CV"
  }
  
  y.ANT <- ANT$Std_RT
  
  ANT_original<-ANT
}

if (DataSet_train == 'HCP' || DataSet_test == 'HCP')
{
  HCP <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/HCP_HOMs_4_ZTransformed.csv", header=FALSE)
  # HCP <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/HCP_HOMs_4_ZTransformed_noABS.csv", header=FALSE)
  # HCP <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/HCP_HOMs_WeightedFC_4_ZTransformed.csv", header=FALSE)
  str(HCP)
  colnames(HCP)[1] <- "individual"
  colnames(HCP)[8:ncol(HCP)] <- c("median RT_0bk (body)","median RT_0bk (face)","median RT_0bk (place)","median RT_0bk (tool)",
                                  "median RT target_0bk (body)","median RT target_0bk (face)", "median RT target_0bk (place)","median RT target_0bk (tool)",
                                  "median RT nontarget_0bk (body)","median RT nontarget_0bk (face)","median RT nontarget_0bk (place)","median RT nontarget_0bk (tool)",
                                  "median RT_2bk (body)","median RT_2bk (face)","median RT_2bk (place)","median RT_2bk (tool)",
                                  "median RT target_2bk (body)","median RT target_2bk (face)", "median RT target_2bk (place)","median RT target_2bk (tool)",
                                  "median RT nontarget_2bk (body)","median RT nontarget_2bk (face)","median RT nontarget_2bk (place)","median RT nontarget_2bk (tool)")
  HCP$medRT<-rowMeans(HCP[,8:31],na.rm=TRUE)
  
  HCP$TC <- rowMeans(cbind(HCP$V2))
  HCP$RTC <- rowMeans(cbind(HCP$V3,HCP$V4,HCP$V5,HCP$V6))
  HCP$RC <- rowMeans(cbind(HCP$V7))
  
  
  runs <-4
  # outliers <-c(which(HCP$medRT %in% boxplot.stats(HCP$medRT,coef = 3)$out))
  behavcol <- grep("medRT", colnames(HCP))
  # outliers <- which(!complete.cases(HCP[,behavcol]))
  outliers <-c(which(HCP$medRT %in% boxplot.stats(HCP$medRT,coef = 3)$out | is.na(HCP$medRT)))
  
  # outliers
  HCP<-HCP[-outliers,]
  
  x.HCP <- HCP[,c(2:7,33:ncol(HCP))]
  y.HCP <- HCP$medRT
  HCP_original<-HCP
  
}


# nback -> gradCPT, ANT, HCP --------------------------------------------------------
if (DataSet_train == 'nback')
{
  nback <-nback_original
  x.nback <- nback[,c(2:11,13:15)]
  y.nback <- nback$dprime
  
  if (DataSet_test == 'gradCPT')
  {
    gradCPT<-gradCPT_original
    gradCPT <- gradCPT[-c(1,3,4,9),]
    x.gradCPT <- gradCPT[,c(2:11,13:15)]
    y.gradCPT <- gradCPT$dprime
    
    if (Analysis == 'HOM')
      # lm.out <- lm(dprime ~ (V2+V3+V6+V11),
      #               data = nback)
      lm.out <- lm(dprime ~ (TC+RC),
                   data = nback)
    if (Analysis == 'Eigen')
      lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5), 
                   data = nback)
    
    lm.TC <- lm(dprime ~ (TC+RC),
                data = nback)
    vif_val <- VIF(lm.TC)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    # lm.out <- lm(dprime ~ (V2+V3+V6+V11+RTC),
    #               data = nback, subset = trainFlag)
    y.glm <- predict(lm.out,gradCPT)
    coeff <- lm.out$coefficients[-1]
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/gradCPT by nback",Analysis,".pdf",sep = "_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.gradCPT)
    correl <- cor.test(y.glm, y.gradCPT, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.gradCPT,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5), cex.axis = 2.5,
         xlab="Observed d'", ylab="Predicted d'", cex.lab = 2.5,
         # xlab="", ylab="", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
    
    
  }
  
  if (DataSet_test == 'ANT')
  {
    ANT <- ANT_original
    ANT <- ANT[-c(1,3,4,13,36),]
    x.ANT <- ANT[,c(2:11,17:19)]
    y.ANT <- ANT$Std_RT
    
    if (Analysis == 'HOM')
      # lm.out <- lm(dprime ~ (V2+V3+V6+V11+RTC),
      # data = nback)
      lm.out <- lm(dprime ~ (TC+RC),
                   data = nback)
    if (Analysis == 'Eigen')
      lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5),
                   data = nback)
    
    
    lm.TC <- lm(dprime ~ (TC+RC),
                data = nback)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,ANT)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/ANT",runs, ANTbehav,"by nback",Analysis,".pdf",sep="_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.ANT)
    correl <- cor.test(y.glm, y.ANT, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.ANT,y.glm,  ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5), cex.axis = 2.5,
         xlab=paste("Observed",behav_label,sep=" "), ylab= "Predicted d'", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
  if (DataSet_test == 'HCP')
  {
    HCP<-HCP_original
    x.HCP <- HCP[,c(2:7,33:ncol(HCP))]
    y.HCP <- HCP$medRT
    
    if (Analysis == 'HOM')
      lm.out <- lm(dprime ~ (TC+RC),
                   data = nback)
    if (Analysis == 'Eigen')
      lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5), 
                   data = nback)
    
    lm.TC <- lm(dprime ~ (TC+RC),
                data = nback)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,HCP)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/HCP by nback",Analysis,".pdf",sep = "_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.HCP)
    correl <- cor.test(y.glm, y.HCP, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.HCP,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5),  cex.axis = 2.5,
         xlab="Observed Median RT", ylab="Predicted d'", cex.lab = 2.5,
         # xlab="", ylab="", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
}


# gradCPT -> nback, ANT, HCP --------------------------------------------------------

if (DataSet_train == 'gradCPT')
{
  
  gradCPT<-gradCPT_original
  x.gradCPT <- gradCPT[,c(2:11,13:15)]
  y.gradCPT <- gradCPT$dprime
  
  if (DataSet_test == 'nback')
  {
    nback <-nback_original
    nback <- nback[-c(1,4,5,9),]
    x.nback <- nback[,c(2:11,13:15)]
    y.nback <- nback$dprime
    
    if (Analysis == 'HOM')
      # lm.out <- lm(dprime ~ (V2+V3+V6+V11),
      #               data = gradCPT)
      lm.out <- lm(dprime ~ (TC+RC),
                   data = gradCPT)
    if (Analysis == 'Eigen')
      lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5),
                   data = gradCPT)
    
    
    lm.TC <- lm(dprime ~ (TC+RC),
                data = gradCPT)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    
    y.glm <- predict(lm.out,nback)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/nback by gradCPT",Analysis,".pdf",sep = "_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.nback)
    correl <- cor.test(y.glm, y.nback, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.nback,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5),  cex.axis = 2.5,
         xlab="Observed d'", ylab="Predicted d'", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
  if (DataSet_test == 'ANT')
  {
    
    ANT <- ANT_original
    ANT <- ANT[-c(1,3,22),]
    x.ANT <- ANT[,c(2:11,17:19)]
    y.ANT <- ANT$Std_RT
    
    if (Analysis == 'HOM')
      # lm.out <- lm(dprime ~ (V2+V3+V6+V11+RTC),
      # data = gradCPT)
      lm.out <- lm(dprime ~ (TC+RC),
                   data = gradCPT)
    if (Analysis == 'Eigen')
      lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5),
                   data = gradCPT)
    
    lm.TC <- lm(dprime ~ (TC+RC),
                data = gradCPT)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,ANT)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/ANT",runs, ANTbehav,"by gradCPT",Analysis,".pdf",sep="_"), height=8, width=8)
    
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.ANT)
    correl <- cor.test(y.glm, y.ANT, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.ANT,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5),  cex.axis = 2.5,
         xlab=paste("Observed",behav_label,sep=" "), ylab="Predicted d'", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
  if (DataSet_test == 'HCP')
  {
    HCP<-HCP_original
    x.HCP <- HCP[,c(2:7,33:ncol(HCP))]
    y.HCP <- HCP$medRT
    
    if (Analysis == 'HOM')
      lm.out <- lm(dprime ~ (TC+RC),
                   data = gradCPT)
    if (Analysis == 'Eigen')
      lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5), 
                   data = nback)
    
    lm.TC <- lm(dprime ~ (TC+RC),
                data = gradCPT)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,HCP)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/HCP by gradCPT",Analysis,".pdf",sep = "_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.HCP)
    correl <- cor.test(y.glm, y.HCP, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.HCP,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5),  cex.axis = 2.5,
         xlab="Observed Median RT", ylab="Predicted d'", cex.lab = 2.5,
         # xlab="", ylab="", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
}



# ANT -> nback, gradCPT, HCP ------------------------------------------------------------

if (DataSet_train == 'ANT')
{
  ANT <- ANT_original
  x.ANT <- ANT[,c(2:11,17:19)]
  y.ANT <- ANT$Std_RT
  
  if (DataSet_test == 'nback')
  {
    nback <-nback_original
    nback <- nback[-c(1,5,6,22,24),]
    x.nback <- nback[,c(2:11,13:15)]
    y.nback <- nback$dprime
    
    if (Analysis == 'HOM')
      # lm.out <- lm(Std_RT ~ (V2+V3+V6+V11+RTC),
      # data = ANT)
      lm.out <- lm(Std_RT ~ (TC+RC),
                   data = ANT)
    if (Analysis == 'Eigen')
      lm.out <- lm(Std_RT ~ (eig1+eig2+eig3+eig4+eig5),
                   data = ANT)
    
    lm.TC <- lm(Std_RT ~ (TC+RC),
                data = ANT)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,nback)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/nback by ANT",runs, ANTbehav,Analysis,".pdf",sep="_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.nback)
    correl <- cor.test(y.glm, y.nback, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.nback,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5),  cex.axis = 2.5,
         xlab="Observed d'", ylab=paste("Observed",behav_label,sep=" "), cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
  if (DataSet_test == 'gradCPT')
  {
    gradCPT<-gradCPT_original
    gradCPT <- gradCPT[-c(4,9,14),]
    x.gradCPT <- gradCPT[,c(2:11,13:15)]
    y.gradCPT <- gradCPT$dprime
    
    if (Analysis == 'HOM')
      # lm.out <- lm(Std_RT ~ (V2+V3+V6+V11+RTC),
      # data = ANT)
      lm.out <- lm(Std_RT ~ (TC+RC),
                   data = ANT)
    if (Analysis == 'Eigen')
      lm.out <- lm(Std_RT ~ (eig1+eig2+eig3+eig4+eig5),
                   data = ANT)
    
    lm.TC <- lm(Std_RT ~ (TC+RC),
                data = ANT)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,gradCPT)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/gradCPT by ANT",runs, ANTbehav,Analysis,".pdf",sep="_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.gradCPT)
    correl <- cor.test(y.glm, y.gradCPT, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.gradCPT,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5),  cex.axis = 2.5,
         xlab="Observed d'", ylab=paste("Observed",behav_label,sep=" "), cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
  if (DataSet_test == 'HCP')
  {
    HCP<-HCP_original
    x.HCP <- HCP[,c(2:7,33:ncol(HCP))]
    y.HCP <- HCP$medRT
    
    if (Analysis == 'HOM')
      lm.out <- lm(Std_RT ~ (TC+RC),
                   data = ANT)
    if (Analysis == 'Eigen')
      lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5), 
                   data = nback)
    
    lm.TC <- lm(Std_RT ~ (TC+RC),
                data = ANT)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,HCP)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/HCP by ANT",runs, ANTbehav,Analysis,".pdf",sep="_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.HCP)
    correl <- cor.test(y.glm, y.HCP, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.HCP,y.glm,  ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5), cex.axis = 2.5,
         xlab="Observed Median RT", ylab=paste("Observed",behav_label,sep=" "), cex.lab = 2.5,
         # xlab="", ylab="", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
}


# HCP -> nback, gradCPT, ANT --------------------------------------------------------
if (DataSet_train == 'HCP')
{
  HCP <-HCP_original
  x.HCP <- HCP[,c(2:7,33:ncol(HCP))]
  y.HCP <- HCP$medRT
  
  if (DataSet_test == 'nback')
  {
    nback<-nback_original
    x.nback <- nback[,c(2:11,13:15)]
    y.nback <- nback$dprime
    
    if (Analysis == 'HOM')
      lm.out <- lm(medRT ~ (TC+RC),
                   data = HCP)
    if (Analysis == 'Eigen')
      lm.out <- lm(medRT ~ (eig1+eig2+eig3+eig4+eig5), 
                   data = HCP)
    
    lm.TC <- lm(medRT ~ (TC+RC),
                data = HCP)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,nback)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/nback by HCP",Analysis,".pdf",sep = "_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.nback)
    correl <- cor.test(y.glm, y.nback, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.nback,y.glm,  ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5), cex.axis = 2.5,
         xlab="Observed d'", ylab="Predicted Median RT", cex.lab = 2.5,
         # xlab="", ylab="", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
  if (DataSet_test == 'gradCPT')
  {
    gradCPT<-gradCPT_original
    x.gradCPT <- gradCPT[,c(2:11,13:15)]
    y.gradCPT <- gradCPT$dprime
    
    if (Analysis == 'HOM')
      lm.out <- lm(medRT ~ (TC+RC),
                   data = HCP)
    if (Analysis == 'Eigen')
      lm.out <- lm(medRT ~ (eig1+eig2+eig3+eig4+eig5), 
                   data = HCP)
    
    lm.TC <- lm(medRT ~ (TC+RC),
                data = HCP)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,gradCPT)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/gradCPT by HCP",Analysis,".pdf",sep = "_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.gradCPT)
    correl <- cor.test(y.glm, y.gradCPT, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.gradCPT,y.glm,  ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5), cex.axis = 2.5,
         xlab="Observed d'", ylab="Predicted Median RT", cex.lab = 2.5,
         # xlab="", ylab="", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
  
  if (DataSet_test == 'ANT')
  {
    ANT <- ANT_original
    ANT <- ANT[-c(1,3,4,13,36),]
    x.ANT <- ANT[,c(2:11,17:19)]
    y.ANT <- ANT$Std_RT
    
    if (Analysis == 'HOM')
      lm.out <- lm(medRT ~ (TC+RC),
                   data = HCP)
    if (Analysis == 'Eigen')
      lm.out <- lm(medRT ~ (eig1+eig2+eig3+eig4+eig5),
                   data = HCP)
    
    lm.TC <- lm(medRT ~ (TC+RC),
                data = HCP)
    anov <- anova(lm.TC,lm.out, test = "Chisq")
    
    y.glm <- predict(lm.out,ANT)
    coeff <- lm.out$coefficients[-1]
    vif_val <- VIF(lm.TC)
    
    pdf(paste("/Users/Mehraveh/Desktop/Final/cross_datasets/ANT",runs, ANTbehav,"by HCP",Analysis,".pdf",sep="_"), height=8, width=8)
    par(pin = c(8,8), mai = rep(1,4))
    reg <- lm(y.glm ~ y.ANT)
    correl <- cor.test(y.glm, y.ANT, method = "pearson", alternative = "two.sided")
    yrange <- (range(y.glm,na.rm = TRUE)[2]-range(y.glm,na.rm = TRUE)[1])
    plot(y.ANT,y.glm, ylim = c(min(y.glm,na.rm=TRUE)-yrange/50,max(y.glm,na.rm=TRUE)+yrange/5),  cex.axis = 2.5,
         xlab=paste("Observed",behav_label,sep=" "), ylab= "Predicted Median RT", cex.lab = 2.5,
         pch=19, cex=2.5)
    options("scipen"=100, digits = 2)
    options(scipen=0,digits = 2)
    legend(x = "topleft", cex = 2.5,
           legend = substitute(list(r == r2, p < p2), list(r2=correl$estimate, p2=format((correl$p.value),scientific=TRUE))))
    abline(reg, col="red", lwd = 7)
    dev.off()
  }
}

normality <- jarque.bera.test(y.glm)$p.value
print(normality)
print(correl$estimate)
print(vif_val)

sum.lm<-summary(lm.out)
print(sum.lm$coefficients[2,4]/2)
print(sum.lm$coefficients[3,4]/2)

writeMat(sprintf(paste("/Users/Mehraveh/Documents/MATLAB/Connectivity/R_results/GLM_coefficients_train_",DataSet_train,"_test_",DataSet_test,"_",runs,".mat",sep = "")), coeff=coeff)
