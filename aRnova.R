library(tibble)
library(dplyr)
library(data.table)
library(ggpubr)
library(stats)

readData <- function(mainDir, subDir){
  fileDir <- paste(mainDir, subDir,sep="")
  spdf <- as.data.frame(read.csv(fileDir, header = T, row.names = 1, stringsAsFactors = FALSE))
  t_spdf <- transpose(spdf)
  colnames(t_spdf) <- rownames(spdf)
  rownames(t_spdf) <- colnames(spdf)
  return(t_spdf)
}

summarizePeaks <- function(spPeak, mainDF){
  #group data by their sample type, requires dplyr for the pipe operation, and perform summary statistics
  #plot distribution of points for visualization
  
  group_by(mainDF, Group) %>%
    summarise(
      count = n(),
      mean = mean(as.vector(spPeak), na.rm = TRUE),
      sd = sd(spPeak, na.rm = TRUE),
      median = median(spPeak, na.rm = TRUE),
      IQR = IQR(spPeak, na.rm = TRUE)
    )
  
  ggboxplot(mainDF, x = "Group", y = as.character(spPeak),
            color = "Group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
            order = c("NC", "SF", "SN"),
            ylab = as.character(spPeak), xlab = "Group")
}

createLinearModel <- function(spPeak, mainDF){
  model <- lm(as.formula(paste(spPeak, " ~ Group",sep="")), data = mainDF)
  res <- resid(model)
  qqnorm(res)
  qqline(res)
  plot(density(res))
}
##################################################
#          MAIN SCRIPT BELOW                     #
##################################################

#read data
spdf <- readData("D:/R scripts/_PPMF-Univariate","/IR Data.csv")

#convert character-formatted numbers to numeric format
spdf[,2:ncol(spdf)] <- lapply(spdf[,2:ncol(spdf)],as.numeric)

#reload because i have no idea why the dataframe is malfunctioning
write.csv(spdf,"D:/R scripts/_PPMF-Univariate/IR Data2.csv")
spdf <- as.data.frame(read.csv("D:/R scripts/_PPMF-Univariate/IR Data2.csv", header = T, row.names = 1, stringsAsFactors = FALSE))

# use for performing summary statistics and checking residuals for normality
# summarizePeak("X499.821364", spdf)
# createLinearModel("X499.821364", spdf)

#initialize vectors for p-values and statistics
p_val_AOV <- c("AOV p-value")
F_val_AOV <- c("AOV F-value")
p_val_KWT <- c("KWT p-value")
Xs_KWT <- c("KWT Chi-squared")

#perform parametric one-way ANOVA and non-parametric kruskal wallis for all peaks
for(i in colnames(spdf)){
  
  if(i == "Group"){
    next
  }
  
  else{
    spdf.aov <- aov(lm(as.formula(paste(i, " ~ Group", sep="")), data = spdf))
    spdf.kwt <- kruskal.test(as.formula(paste(i, " ~ Group", sep="")), data = spdf)
    
    # flig <- fligner.test(as.formula(paste(i, " ~ Group", sep="")), data = spdf)
    # if(flig$p.value < 0.05){
    #   print(flig)
    # }
    
    summary(spdf.aov)
    p_val_AOV <- append(p_val_AOV, summary(spdf.aov)[[1]][["Pr(>F)"]][1])
    F_val_AOV <- append(F_val_AOV, summary(spdf.aov)[[1]][["F value"]][1])
    p_val_KWT <- append(p_val_KWT, spdf.kwt$p.value)
    Xs_KWT <- append(Xs_KWT, spdf.kwt$statistic)
  }
}

#summarize and make a result table, create a csv output
spdf_res <- do.call("rbind", list(spdf, p_val_AOV, F_val_AOV, p_val_KWT, Xs_KWT))
spdf_res <- t(spdf_res)
write.csv(spdf_res, "D:/R scripts/_PPMF-Univariate/IR Data2.csv")

