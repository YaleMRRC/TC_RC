rm(list = ls())
library(class)
library(FNN)
library(gtools)
library(R.matlab)
library(multcomp)
library(MVN)
library(moments)
library(tseries)

normality <- matrix(0,4,1)
# Analysis <- readline(cat("What type of features? ('HOM','Eigen')\n"))
Analysis <- 'HOM'
DataSet <- readline(cat("What data set? ('nback', 'gradCPT', 'ANT', 'both (nback & gradCPT)', 'HCP' )\n"))

options(digits=2)
dir.create("/Users/Mehraveh/Desktop/Final/",showWarnings = FALSE)
dir.create("/Users/Mehraveh/Desktop/Final/within_datasets/",showWarnings = FALSE)

if (DataSet == 'nback' || DataSet == 'both')
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



if (DataSet == 'gradCPT' || DataSet == 'both')
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



if (DataSet == 'ANT')
{
  # runs <- readline(cat("What size HOM? ('5', '7')\n"))
  runs <- 7
  runs <- as.integer(runs)
  
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

if (DataSet == 'HCP')
{
  HCP <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/HCP_HOMs_4_ZTransformed.csv", header=FALSE)
  # HCP <- read.csv("/Users/Mehraveh/Documents/MATLAB/Connectivity/NetworkMeasures/ZTransformed/HCP_HOMs_WeightedFC_4_ZTransformed.csv", header=FALSE)
  str(HCP)
  colnames(HCP)[1] <- "individual"
  colnames(HCP)[8:ncol(HCP)] <- c("median RT_0bk (body)","median RT_0bk (face)","median RT_0bk (place)","median RT_0bk (tool)",
                                  "median RT target_0bk (body)","median RT target_0bk (face)", "median RT target_0bk (place)","median RT target_0bk (tool)",
                                  "median RT nontarget_0bk (body)","median RT nontarget_0bk (face)","median RT nontarget_0bk (place)","median RT nontarget_0bk (tool)",
                                  "median RT_2bk (body)","median RT_2bk (face)","median RT_2bk (place)","median RT_2bk (tool)",
                                  "median RT target_2bk (body)","median RT target_2bk (face)", "median RT target_2bk (place)","median RT target_2bk (tool)",
                                  "median RT nontarget_2bk (body)","median RT nontarget_2bk (face)","median RT nontarget_2bk (place)","median RT nontarget_2bk (tool)")
  HCP$medRT<-rowMeans(HCP[,8:ncol(HCP)],na.rm=TRUE)
  
  HCP$TC <- rowMeans(cbind(HCP$V2))
  HCP$RTC <- rowMeans(cbind(HCP$V3,HCP$V4,HCP$V5,HCP$V6))
  HCP$RC <- rowMeans(cbind(HCP$V7))
  
  
  runs <-4
  outliers <-c(which(HCP$medRT %in% boxplot.stats(HCP$medRT,coef = 3)$out | is.na(HCP$medRT)))
  outliers
  HCP<-HCP[-outliers,]
  x.HCP <- HCP[,c(2:7,33:ncol(HCP))]
  y.HCP <- HCP$medRT
  
  HCP_original<-HCP
  
}


# nback -------------------------------------------------------------------
if (DataSet == 'nback')
{
  nback<-nback_original
  dprime_original <- nback$dprime
  x <- x.nback
  l = nrow(x)
  perm_bool <- menu(c("Yes", "No"), title="Do you want permutation test?")
  
  if (perm_bool==1){
    P <- 100
  } else {
    P <- 1
  }
  
  # pvalue = 9.9e-03
  nn=0
  for (perm in 1:P)
  {
    print(perm)
    if (perm == 1){
      nback$dprime <- dprime_original
    }else{
      nback$dprime <- permute(dprime_original)
    }
    
    y <- nback$dprime
    
    anov <- matrix(0,l,1)
    y.lm <- matrix(0,l,1)
    coeff <- matrix(0,l,2)
    for (i in 1:l)
    {
      testFlag = logical(length = l)
      testFlag[i] = TRUE
      trainFlag <- !testFlag
      xtest <- x[!trainFlag,]
      ytest <- y[!trainFlag]
      
      if (Analysis == 'HOM')
        # lm.out <- lm(dprime ~ (V2+V3+V6+V11),
        #               data = nback, subset = trainFlag)
        lm.out <- glm(dprime ~ (TC+RC),
                      data = nback, subset = trainFlag, family=gaussian())
      if (Analysis == 'Eigen')
        lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5),
                     data = nback, subset = trainFlag)
      
      
      lm.TC <- glm(dprime ~ (TC+RC),
                   data = nback, subset = trainFlag)
      tmp <- anova(lm.TC,lm.out, test = "F")
      anov[i] <- tmp$`Pr(>F)`[2]
      
      y.lm [i] <- predict(lm.out,xtest)
      coeff[i,] <- lm.out$coefficients[-1]
      
      
    }
    reg <- lm((y.lm) ~ (y))
    correl <- cor.test(y.lm, y, method = "pearson", alternative = "two.sided")
    if (perm == 1){
      r.original = correl$estimate
      y.original = y
      y.lm.original = y.lm
      print(r.original)
    } else {
      r = correl$estimate
      if (r >= r.original)
        nn=nn+1
    }
  }
  
  pvalue = (nn+1)/(100+1)
  
  pdf(paste("/Users/Mehraveh/Desktop/Final/within_datasets/nback",Analysis,".pdf", sep = "_"), height=8, width=8)
  par(pin = c(8,8), mai = rep(1,4))
  reg <- lm(y.lm.original ~ y.original)
  yrange <- (range(y.lm.original,na.rm = TRUE)[2]-range(y.lm.original,na.rm = TRUE)[1])
  plot(y.original,y.lm.original,ylim = c(min(y.lm.original,na.rm=TRUE)-yrange/50,max(y.lm.original,na.rm=TRUE)+yrange/5), cex.axis = 2.5,
       xlab="Observed d'", ylab="Predicted d'", cex.lab = 2.5,
       pch=19, cex=2.5)
  options("scipen"=-100, digits = 2)
  options(scipen=0, digits = 2)
  legend(x = "topleft", cex = 2.5,
         legend = substitute(list(r == r2, p < p2), list(r2=r.original, p2=format(pvalue,scientific=TRUE))))
  abline(reg, col="red", lwd = 7)
  dev.off()
  
  # ANOVA and Dunnett's test
  group <- c(rep("TC1",l), rep("TC2",l), rep("TC3",l), rep("RC",l))
  value <- c(nback$V2,nback$V3,nback$V6,nback$RC)
  data <- data.frame(value,group)
  aov1 <- aov(value ~ group, data=data)
  set.seed(1)
  results<-summary(glht(aov1, linfct=mcp(group="Dunnett"),alternative = c("two.sided")))
  pval=results$test$pvalues[1:sum(1:(runs-3))]
  
  
  # simple t-test to comprea TC_avg and RC
  y <- c(nback$TC,nback$RC)
  group <- c(rep(1, times=c(length(nback$TC))),rep(2, times=c(length(nback$RC))))
  
  h<-pairwise.t.test(y, group, alternative = "less")
  h<-wilcox.test(y,group,alternative = "less")
  pval <- h$p.value
  #pval=pval[1,1]
  # normality <- mardiaTest(data.frame(nback$TC,nback$RC), cov = TRUE, qqplot = TRUE)
  normality[1] <- jarque.bera.test(nback$TC)$p.value
  normality[2] <- jarque.bera.test(nback$RC)$p.value
  normality[3] <- jarque.bera.test(y.original)$p.value
  normality[4] <- jarque.bera.test(as.vector(y.lm.original))$p.value
  
  
}

# gradCPT -----------------------------------------------------------------

if (DataSet == 'gradCPT')
{
  gradCPT<-gradCPT_original
  dprime_original <- gradCPT$dprime
  x <- x.gradCPT
  l = nrow(x)
  
  perm_bool <- menu(c("Yes", "No"), title="Do you want permutation test?")
  
  
  if (perm_bool==1){
    P <- 100
  } else {
    P <- 1
  }
  
  # pvalue = 9.9e-03
  nn=0
  for (perm in 1:P)
  {
    print(perm)
    if (perm == 1){
      gradCPT$dprime <- dprime_original
    }else{
      gradCPT$dprime <- permute(dprime_original)
    }
    
    y <- gradCPT$dprime
    coeff <- matrix(0,l,2)
    anov <- matrix(0,l,1)
    y.lm <- matrix(0,l,1)
    for (i in 1:l)
    {
      testFlag = logical(length = l)
      testFlag[i] = TRUE
      trainFlag <- !testFlag
      xtest <- x[!trainFlag,]
      ytest <- y[!trainFlag]
      
      if (Analysis == 'HOM')
        # lm.out <- lm(dprime ~ (V2+V3+V6+V11),
        #             data = gradCPT, subset = trainFlag)
        lm.out <- lm(dprime ~ (TC+RC),
                     data = gradCPT, subset = trainFlag, family=gaussian())
      
      
      if (Analysis == 'Eigen')
        lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5),
                     data = gradCPT, subset = trainFlag)
      
      lm.TC <- lm(dprime ~ (TC+RC),
                  data = gradCPT, subset = trainFlag)
      tmp <- anova(lm.TC,lm.out, test = "F")
      anov[i] <- tmp$`Pr(>F)`[2]
      
      y.lm [i] <- predict(lm.out,xtest)
      coeff[i,] <- lm.out$coefficients[-1]
    }
    
    reg <- lm((y.lm) ~ (y))
    correl <- cor.test(y.lm, y, method = "pearson", alternative = "two.sided")
    if (perm == 1){
      r.original = correl$estimate
      y.original = y
      y.lm.original = y.lm
      print(r.original)
    } else {
      r = correl$estimate
      if (r >= r.original)
        nn=nn+1
    }
  }
  
  pvalue = (nn+1)/(100+1)
  
  
  pdf(paste("/Users/Mehraveh/Desktop/Final/within_datasets/gradCPT",Analysis,".pdf",sep = "_"), height=8, width=8)
  par(pin = c(8,8), mai = rep(1,4))
  reg <- lm(y.lm.original ~ y.original)
  yrange <- (range(y.lm.original,na.rm = TRUE)[2]-range(y.lm.original,na.rm = TRUE)[1])
  plot(y.original,y.lm.original,ylim = c(min(y.lm.original,na.rm = TRUE)-yrange/50,max(y.lm.original,na.rm = TRUE)+yrange/5), cex.axis = 2.5, 
       xlab="Observed d'", ylab="Predicted d'", cex.lab = 2.5,
       # xlab="", ylab="", cex.lab = 2.5, 
       pch=19, cex=2.5)
  options("scipen"=100, digits = 2)
  options(scipen=0, digits = 2)
  legend(x = "topleft", cex = 2.5,
         legend = substitute(list(r == r2, p < p2), list(r2=r.original, p2=format(pvalue,scientific=TRUE))))
  abline(reg, col="red", lwd=7)
  dev.off()
  
  # ANOVA and Dunnett's test
  group <- c(rep("TC1",l), rep("TC2",l), rep("TC3",l), rep("RC",l))
  value <- c(gradCPT$V2,gradCPT$V3,gradCPT$V6,gradCPT$RC)
  data <- data.frame(value,group)
  aov1 <- aov(value ~ group, data=data)
  set.seed(1)
  results<-summary(glht(aov1, linfct=mcp(group="Dunnett"),alternative = c("two.sided")))
  pval=results$test$pvalues[1:sum(1:(runs-3))]
  
  # simple t-test to comprea TC_avg and RC
  TC_RC <- c(gradCPT$TC,gradCPT$RC)
  group <- c(rep(1, times=c(length(gradCPT$TC))),rep(2, times=c(length(gradCPT$RC))))
  
  h<-pairwise.t.test(TC_RC, group, alternative = "less")
  h<-wilcox.test(TC_RC,group,alternative = "less")
  pval <- h$p.value
  #pval=pval[1,1]
  # normality <- mardiaTest(data.frame(gradCPT$TC,gradCPT$RC), cov = TRUE, qqplot = TRUE)
  normality[1] <- jarque.bera.test(gradCPT$TC)$p.value
  normality[2] <- jarque.bera.test(gradCPT$RC)$p.value
  normality[3] <- jarque.bera.test(y.original)$p.value
  normality[4] <- jarque.bera.test(as.vector(y.lm.original))$p.value
  
  
}


# ANT ---------------------------------------------------------------------
if (DataSet == 'ANT')
{
  ANT<-ANT_original
  dprime_original <- ANT$Std_RT
  x <- x.ANT
  l = nrow(x)
  
  perm_bool <- menu(c("Yes", "No"), title="Do you want permutation test?")
  
  if (perm_bool==1){
    P <- 100
  } else {
    P <- 1
  }
  
  
  nn=0
  for (perm in 1:P)
  {
    print(perm)
    if (perm == 1){
      ANT$Std_RT <- dprime_original
    }else{
      ANT$Std_RT <- permute(dprime_original)
    }
    y <- ANT$Std_RT
    
    
    coeff <- matrix(0,l,2)
    anov <- matrix(0,l,1)
    y.lm <- matrix(0,l,1)
    for (i in 1:l)
    {
      testFlag = logical(length = l)
      testFlag[i] = TRUE
      trainFlag <- !testFlag
      xtest <- x[!trainFlag,]
      ytest <- y[!trainFlag]
      
      if (runs == 5)
      {
        
        if (Analysis == 'HOM')
        {
          # lm.out <- lm(Std_RT ~ (V2+V3+V6+V11),
          # data = ANT, subset = trainFlag)
          lm.out <- lm(Std_RT ~ (TC+RC),
                       data = ANT, subset = trainFlag, family=gaussian)
          
        }
        if (Analysis == 'Eigen')
          lm.out <- lm(Std_RT ~ (eig1+eig2+eig3+eig4+eig5),
                       data = ANT, subset = trainFlag)
        
      }
      
      if (runs == 7)
      {
        if (Analysis == 'HOM')
          # lm.out <- lm(Std_RT ~ (V2+V3+V4+V5+V8+V9+V10+V13+V14+V17+V22), #+RTC),
          #               data = ANT, subset = trainFlag)
          lm.out <- lm(Std_RT ~ (TC+RC),
                       data = ANT, subset = trainFlag)
        
        if (Analysis == 'Eigen')
          lm.out <- lm(Std_RT ~ (eig1+eig2+eig3+eig4+eig5+eig6+eig7),
                       data = ANT, subset = trainFlag)
        
      }
      
      lm.TC <- lm(Std_RT ~ (TC+RC),
                  data = ANT, subset = trainFlag)
      tmp <- anova(lm.TC,lm.out, test = "F")
      anov[i] <- tmp$`Pr(>F)`[2]
      
      y.lm [i] <- predict(lm.out,xtest)
      coeff[i,] <- lm.out$coefficients[-1]
    }
    
    reg <- lm((y.lm) ~ (y))
    correl <- cor.test(y.lm, y, method = "pearson", alternative = "two.sided")
    if (perm == 1){
      r.original = correl$estimate
      y.original = y
      y.lm.original = y.lm
      print(r.original)
    } else {
      r = correl$estimate
      if (r > r.original)
        nn=nn+1
    }
  }
  
  pvalue = (nn+1)/(100+1)
  
  
  pdf(paste("/Users/Mehraveh/Desktop/Final/within_datasets/ANT",runs,ANTbehav,Analysis,".pdf",sep="_"), height=8, width=8)
  par(pin = c(8,8), mai = rep(1,4))
  reg <- lm(y.lm.original ~ y.original)
  yrange <- (range(y.lm.original,na.rm = TRUE)[2]-range(y.lm.original,na.rm = TRUE)[1])
  plot(y.original,y.lm.original,ylim = c(min(y.lm.original,na.rm = TRUE)-yrange/50,max(y.lm.original,na.rm = TRUE)+yrange/5),cex.axis = 2.5,
       xlab=paste("Observed",behav_label,sep=" "), ylab=paste("Predicted",behav_label,sep=" ") , cex.lab = 2.5,
       pch=19, cex=2.5)
  options("scipen"=-100, digits = 2)
  options(scipen=0, digits = 2)
  legend(x = "topleft", cex = 2.5,
         legend = substitute(list(r == r2, p < p2), list(r2=r.original, p2=format(pvalue,scientific=TRUE))))
  abline(reg, col="red", lwd=7)
  dev.off()
  
  # ANOVA and Dunnett's test
  group <- c(rep("TC1",l), rep("TC2",l), rep("TC3",l),rep("TC4",l),rep("TC5",l),rep("TC6",l),rep("TC7",l),rep("TC8",l), rep("TC9",l),rep("TC10",l),rep("RC",l))
  value <- c(ANT$V2,ANT$V3,ANT$V4,ANT$V5,ANT$V8,ANT$V9,ANT$V10,ANT$V13,ANT$V14,ANT$V17,ANT$RC)
  data <- data.frame(value,group)
  aov1 <- aov(value ~ group, data=data)
  
  results<-summary(glht(aov1, linfct=mcp(group="Dunnett"),alternative = c("two.sided")))
  results
  pval=results$test$pvalues[1:sum(1:(runs-3))]
  
  # simple t-test to comprea TC_avg and RC
  y <- c(ANT$TC,ANT$RC)
  group <- c(rep(1, times=c(length(ANT$TC))),rep(2, times=c(length(ANT$RC))))
  
  h<-pairwise.t.test(y, group, alternative = "less")
  h<-wilcox.test(y,group,alternative = "less")
  pval <- h$p.value
  #pval=pval[1,1]
  # normality <- mardiaTest(data.frame(ANT$TC,ANT$RC), cov = TRUE, qqplot = TRUE)
  normality[1] <- jarque.bera.test(ANT$TC)$p.value
  normality[2] <- jarque.bera.test(ANT$RC)$p.value
  normality[3] <- jarque.bera.test(y.original)$p.value
  normality[4] <- jarque.bera.test(as.vector(y.lm.original))$p.value
  
}





# both (nback + gradcPT) ----------------------------------------------------
if (DataSet == 'both')
{
  
  gradCPT <- gradCPT_original
  gradCPT <- gradCPT[-c(1,3),]
  x.gradCPT <- gradCPT[,c(2:11,13:ncol(gradCPT))]
  y.gradCPT <- gradCPT$dprime
  
  nback<-nback_original
  nback <- nback[-c(1,5),]
  x.nback <- nback[,c(2:11,13:ncol(nback))]
  y.nback <- nback$dprime
  
  
  total <- rbind(gradCPT,nback)
  x <- rbind(x.gradCPT, x.nback)
  y <- total$dprime
  dprime_original <- total$dprime
  l = nrow(total)
  
  
  perm_bool <- menu(c("Yes", "No"), title="Do you want permutation test?")
  
  
  if (perm_bool==1){
    P <- 100
  } else {
    P <- 1
  }
  
  # pvalue = 9.9e-03
  nn=0
  for (perm in 1:P)
  {
    print(perm)
    if (perm == 1){
      total$dprime <- dprime_original
    }else{
      total$dprime <- permute(dprime_original)
    }
    y <- total$dprime
    
    coeff <- matrix(0,l,2)
    anov <- matrix(0,l,1)
    y.lm <- matrix(0,l,1)
    for (i in 1:l)
    {
      
      testFlag = logical(length = l)
      testFlag[i] = TRUE
      trainFlag <- !testFlag
      xtest <- x[!trainFlag,]
      ytest <- y[!trainFlag]
      
      if (Analysis == 'HOM')
        # lm.out <- lm(dprime ~ (V2+V3+V6+V11),
        #               data = total, subset = trainFlag)
        lm.out <- lm(dprime ~ (TC+RC),
                     data = total, subset = trainFlag)
      if (Analysis == 'Eigen')
        lm.out <- lm(dprime ~ (eig1+eig2+eig3+eig4+eig5),
                     data = total, subset = trainFlag)
      
      y.lm[i,] <- predict(lm.out, xtest)
      coeff[i,] <- lm.out$coefficients[-1]
    }
    reg <- lm((y.lm) ~ (y))
    correl <- cor.test(y.lm, y, method = "pearson", alternative = "two.sided")
    if (perm == 1){
      r.original = correl$estimate
      y.original = y
      y.lm.original = y.lm
      print(r.original)
    } else {
      r = correl$estimate
      if (r>= r.original)
        nn=nn+1
    }
  }
  
  pvalue = (1+nn)/(100+1)
  
  
  pdf(paste("/Users/Mehraveh/Desktop/Final/within_datasets/both_leaveoneout",Analysis,".pdf",sep = "_"), height=8, width=8)
  par(pin = c(8,8), mai = rep(1,4))
  reg <- lm(y.lm.original ~ y.original)
  yrange <- (range(y.lm.original,na.rm = TRUE)[2]-range(y.lm.original,na.rm = TRUE)[1])
  plot(y.original[1:14],y.lm.original[1:14],ylim = c(min(y.lm.original,na.rm = TRUE)-yrange/50,max(y.lm.original,na.rm = TRUE)+yrange/5), cex.axis = 2.5,
       xlab="Observed d'", ylab="Predicted d'", cex.lab = 2.5,
       pch=19, cex=2.5)
  points(y.original[15:39],y.lm.original[15:39],pch=2, cex = 2.5)
  options("scipen"=100, digits = 2)
  options(scipen=0, digits = 2)
  legend(x = "topleft", cex = 2.5,
         legend = substitute(list(r == r2, p < p2), list(r2=r.original,p2=format(pvalue,scientific=TRUE))))
  legend(x = "bottomright", cex = 2.5,
         c("gradCPT","nback"),pch=c(19,2))
  abline(reg, col="red", lwd = 7)
  dev.off()
  
}



# HCP -------------------------------------------------------------------
if (DataSet == 'HCP')
{
  HCP<-HCP_original
  medRT_original <- HCP$medRT
  x <- x.HCP
  l = nrow(x)
  perm_bool <- menu(c("Yes", "No"), title="Do you want permutation test?")
  
  if (perm_bool==1){
    P <- 100
  } else {
    P <- 1
  }
  
  # pvalue = 9.9e-03
  nn=0
  for (perm in 1:P)
  {
    print(perm)
    if (perm == 1){
      HCP$medRT <- medRT_original
    }else{
      HCP$medRT <- permute(medRT_original)
    }
    
    y <- HCP$medRT
    
    coeff <- matrix(0,l,2)
    anov <- matrix(0,l,1)
    y.lm <- matrix(0,l,1)
    for (i in 1:l)
    {
      testFlag = logical(length = l)
      testFlag[i] = TRUE
      trainFlag <- !testFlag
      xtest <- x[!trainFlag,]
      ytest <- y[!trainFlag]
      
      if (Analysis == 'HOM')
        lm.out <- lm(medRT ~ (TC+RC),
                     data = HCP, subset = trainFlag, family=gaussian)
      if (Analysis == 'Eigen')
        lm.out <- lm(medRT ~ (eig1+eig2+eig3+eig4+eig5),
                     data = HCP, subset = trainFlag)
      
      lm.TC <- lm(medRT ~ (TC+RC),
                  data = HCP, subset = trainFlag)
      tmp <- anova(lm.TC,lm.out, test = "F")
      anov[i] <- tmp$`Pr(>F)`[2]
      
      y.lm [i] <- predict(lm.out,xtest)
      coeff[i,] <- lm.out$coefficients[-1]
    }
    reg <- lm((y.lm) ~ (y))
    correl <- cor.test(y.lm, y, method = "pearson", alternative = "two.sided")
    if (perm == 1){
      r.original = correl$estimate
      y.original = y
      y.lm.original = y.lm
      print(r.original)
    } else {
      r = correl$estimate
      if (r >= r.original)
        nn=nn+1
    }
  }
  
  pvalue = (nn+1)/(100+1)
  
  pdf(paste("/Users/Mehraveh/Desktop/Final/within_datasets/HCP",Analysis,".pdf", sep = "_"), height=8, width=8)
  par(pin = c(8,8), mai = rep(1,4))
  reg <- lm(y.lm.original ~ y.original)
  yrange <- (range(y.lm.original,na.rm = TRUE)[2]-range(y.lm.original,na.rm = TRUE)[1])
  plot(y.original,y.lm.original,ylim = c(min(y.lm.original,na.rm = TRUE)-yrange/50,max(y.lm.original,na.rm = TRUE)+yrange/5), cex.axis = 2.5,
       xlab="Observed Median RT", ylab="Predicted Median RT", cex.lab = 2.5,
       pch=19, cex=2.5)
  options("scipen"=-100, digits = 2)
  options(scipen=0, digits = 2)
  legend(x = "topleft", cex = 2.5,
         legend = substitute(list(r == r2, p < p2), list(r2=r.original, p2=format(pvalue,scientific=TRUE))))
  abline(reg, col="red", lwd = 7)
  dev.off()
  
  # ANOVA and Dunnett's test
  # group <- c(rep("TC",l),rep("RC",l))
  # value <- c(HCP$TC,HCP$RC)
  # data <- data.frame(value,group)
  # aov1 <- aov(value ~ group, data=data)
  # set.seed(1)
  # results<-summary(glht(aov1, linfct=mcp(group="Dunnett"),alternative = c("two.sided")))
  # pval=results$test$pvalues[1:sum(1:(runs-3))]
  
  # simple t-test to comprea TC_avg and RC
  y <- c(HCP$TC,HCP$RC)
  group <- c(rep(1, times=c(length(HCP$TC))),rep(2, times=c(length(HCP$RC))))
  
  h<-pairwise.t.test(y, group, alternative = "less")
  h<-wilcox.test(y,group,alternative = "less")
  pval <- h$p.value
  #pval=pval[1,1]
  # normality <- mardiaTest(data.frame(HCP$TC,HCP$RC), cov = TRUE, qqplot = TRUE)
  normality[1] <- jarque.bera.test(HCP$TC)$p.value
  normality[2] <- jarque.bera.test(HCP$RC)$p.value
  normality[3] <- jarque.bera.test(y.original)$p.value
  normality[4] <- jarque.bera.test(as.vector(y.lm.original))$p.value
  
  
}
# print(normality)
writeMat(sprintf(paste("/Users/Mehraveh/Documents/MATLAB/Connectivity/R_results/ANOVA_Dunnett_Pvals_",DataSet,"_",runs,".mat",sep = "")), pval=pval)  
if (perm_bool==0){
  writeMat(sprintf(paste("/Users/Mehraveh/Documents/MATLAB/Connectivity/R_results/GLM_coefficients_",DataSet,"_",runs,".mat",sep = "")), coeff=coeff)  
}
