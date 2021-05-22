# ### auxiliary functions
# for (i in 1){
#   path <- '/Users/yyq/Library/Mobile Documents/com~apple~CloudDocs/ipad_mac/github/[Yuan2015] BitPhylogeny/我的 test/Summary/'
#   source(paste0(path, 'auxiliary_func.R'), echo=F)
# }

### Read the mnist features extracted from AE
fin = "/Users/yyq/Library/Mobile Documents/com~apple~CloudDocs/ipad_mac/github/[Yuan2015] BitPhylogeny/我的 test/mnist"

mnistDataAE <- read.csv(file = paste(fin, "mnistFeaturesAE5D1.csv", sep = "/"), header = F)
mnistFeaturesAE <- mnistDataAE[,-1]
mnistLabelAE <- mnistDataAE[,1]
# numOfData <- nrow(mnistDataAE)

k = 0
for (i in 0:9) {
  numOfEachDigit <- sum(mnistLabelAE == i)
  cat("Indices of", i, "is:", 
      paste(k+1, k+numOfEachDigit,sep="-"), 
      "(", numOfEachDigit, ",", round(numOfEachDigit/numOfData * 100,1), "% )\n")
  k <- k + numOfEachDigit
}

### miniTrain set: subset obs for each digit
set.seed(fixSeed)

numOfEachDigit <- 100
SelMiniData <- function(i) {
  indices <- sample(which(mnistLabelAE == i), numOfEachDigit)
  mnistDataAE[indices,]
}

miniData <- Reduce(rbind, Map(SelMiniData, 0:9))
mini_X <- miniData[,-1]
mini_y <- miniData[,1]

### 确定用于训练的 testData
# mini-Data
testData <- as.matrix(mini_X)
dim(testData)
numOfData <- nrow(testData)
colnames(testData) <- paste0("V", 1:dim(testData)[2])

p <- ggplot(data = data.frame(x = mini_X[,1], y = mini_X[,2], digits = mini_y), 
            aes(x=x, y=y, group=digits)) +
  geom_point(aes(color = factor(digits)))
p

plot(testData)

# full-Data
testData <- as.matrix(mnistFeaturesAE)
dim(testData)
numOfData <- nrow(testData)
colnames(testData) <- paste0("V", 1:dim(testData)[2])

p <- ggplot(data = data.frame(x = mnistFeaturesAE[,1], y = mnistFeaturesAE[,2], digits = mnistLabelAE), 
            aes(x=x, y=y, group=digits)) +
  geom_point(aes(color = factor(digits)))
p

plot(testData)

