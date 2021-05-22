### DataSet

ImportDataset(fixSeed) {

####### 数据集：seven class (balanced) ########
set.seed(fixSeed)
m <- 700
dims <- 2
testData <- rbind(rmvnorm(m/7, mean = rep(1.0, dims), sigma = diag(0.08^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(1.8, dims), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(m/7, mean = c(2.0, 2.2), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(m/7, mean = c(2.2, 2.0), sigma = diag(0.03^2, dims, dims)),
                  
                  rmvnorm(m/7, mean = rep(0.5, dims), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.5, 0.3), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.3, 0.5), sigma = diag(0.01^2, dims, dims))
)
plot(testData[,1], testData[,2])

#-------------------------------------------------------------------------------
####### 数据集：seven class (balanced & children 方差更小) ########
set.seed(fixSeed)
m <- 700
dims <- 2
testData <- rbind(rmvnorm(m/7, mean = rep(1.0, dims), sigma = diag(0.10^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(1.5, dims), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(m/7, mean = c(1.6, 1.8), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(1.8, 1.6), sigma = diag(0.01^2, dims, dims)),
                  
                  rmvnorm(m/7, mean = rep(0.5, dims), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.5, 0.3), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.3, 0.5), sigma = diag(0.01^2, dims, dims))
)
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集：seven class (balanced & 距离更远) ########
set.seed(fixSeed)
m <- 700
dims <- 2
testData <- rbind(rmvnorm(m/7, mean = rep(2.0, dims), sigma = diag(0.10^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(3.5, dims), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(m/7, mean = c(3.8, 4.0), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(4.0, 3.8), sigma = diag(0.01^2, dims, dims)),
                  
                  rmvnorm(m/7, mean = rep(0.7, dims), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.5, 0.3), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.3, 0.5), sigma = diag(0.01^2, dims, dims))
)
plot(testData[,1], testData[,2])



#-------------------------------------------------------------------------------
####### 数据集：seven class (balanced & 距离更远 & 方差更大) ########
set.seed(fixSeed)
m <- 700
dims <- 2
testData <- rbind(rmvnorm(m/7, mean = rep(2.0, dims), sigma = diag(0.2^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(3.5, dims), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = c(3.8, 4.0), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(4.0, 3.8), sigma = diag(0.01^2, dims, dims)),
                  
                  rmvnorm(m/7, mean = rep(0.7, dims), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.5, 0.3), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(m/7, mean = c(0.3, 0.5), sigma = diag(0.01^2, dims, dims))
)
plot(testData[,1], testData[,2])




#-------------------------------------------------------------------------------
####### 数据集：seven class (Unbalanced) ########
set.seed(fixSeed)
dims <- 2
testData <- rbind(rmvnorm(200, mean = rep(1.0, dims), sigma = diag(0.10^2, dims, dims)),
                  rmvnorm(100, mean = rep(1.5, dims), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(50, mean = c(1.6, 1.8), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(50, mean = c(1.8, 1.6), sigma = diag(0.01^2, dims, dims)),
                  
                  rmvnorm(100, mean = rep(0.5, dims), sigma = diag(0.03^2, dims, dims)),
                  rmvnorm(50, mean = c(0.5, 0.3), sigma = diag(0.01^2, dims, dims)),
                  rmvnorm(50, mean = c(0.3, 0.5), sigma = diag(0.01^2, dims, dims))
)
m <- nrow(testData)
plot(testData[,1], testData[,2])




#-------------------------------------------------------------------------------
####### 数据集 1：两个clusters low-dim ########
set.seed(fixSeed)
m <- 10
dims <- 2
testData <- rbind(matrix(0.05*rnorm(m)+0.4,m/2,dims),
                  matrix(0.01*rnorm(m)+0.6,m/2,dims))
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集 2：两个clusters，high-dim ########
set.seed(fixSeed)
m <- 2000  # 60
dims <- 200
testData <- rbind(matrix(0.05*rnorm(m/2 * dims)+0.5,m/2,dims),
                  matrix(0.01*rnorm(m/2 * dims)+0.0,m/2,dims))
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集 5：两个clusters，medium-dim ########
set.seed(fixSeed)
m <- 200
dims <- 5
testData <- rbind(matrix(0.051*rnorm(m/2 * dims)+0.5,m/2,dims),
                  matrix( 0.01*rnorm(m/2 * dims)+0.0,m/2,dims))
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集 6：两个clusters，multi-normal, medium-dim ########
set.seed(fixSeed)
m <- 30
dims <- 20
mu1 <- rep(2.0, dims)
mu2 <- rep(0.0, dims)
sig1 <- diag(0.05, dims)
diag(sig1[-1,]) <- 0.01
diag(sig1[,-1]) <- 0.01
sig2 <- sig1 * 0.5

testData <- rbind(MASS::mvrnorm(m/2, mu = mu1, Sigma = sig1),
                  MASS::mvrnorm(m/2, mu = mu2, Sigma = sig2))
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集 3：三个clusters low-dim  ########
set.seed(fixSeed)
m <- 300
dims <- 2
testData <- rbind(matrix(0.050*rnorm(m/3 * dims)+2.0,m/3,dims),
                  matrix(0.010*rnorm(m/3 * dims)+3.0,m/3,dims),
                  matrix(0.025*rnorm(m/3 * dims)+1.0,m/3,dims))
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集 4：三个clusters，high-dim ########
# 描述：三个class，有两个距离很近
set.seed(fixSeed)
m <- 3000
dims <- 50
testData <- rbind(matrix(0.05*rnorm(m/3 * dims)+2.0,m/3,dims),
                  matrix(0.01*rnorm(m/3 * dims)+2.2,m/3,dims),
                  matrix(0.02*rnorm(m/3 * dims)+1.0,m/3,dims))
plot(testData[,1], testData[,2])
dim(testData)

#-------------------------------------------------------------------------------
####### 数据集 9：三个clusters，high-dim ########
# 描述：三个class，有两个距离很近，
#       且离得近的两个cluster方差也相近
set.seed(fixSeed)
m <- 3000
dims <- 50
testData <- rbind(matrix(0.05*rnorm(m/3 * dims)+1.0,m/3,dims),
                  matrix(0.01*rnorm(m/3 * dims)+2.2,m/3,dims),
                  matrix(0.02*rnorm(m/3 * dims)+2.0,m/3,dims))
plot(testData[,1], testData[,2])
dim(testData)


#-------------------------------------------------------------------------------
####### 数据集 7：两个cluster，椭圆 ########
set.seed(fixSeed)
m <- 200
dims <- 2
sigma <- matrix(c(.01,.005,.005,.01),2,2)
dat1 <- rmvnorm(m/2, mean = rep(0,dims), sigma = sigma)
dat2 <- rmvnorm(m/2, mean = rep(.7,dims), sigma = sigma)
testData <- rbind(dat1, dat2)
plot(testData[,1], testData[,2])
dim(testData)


#-------------------------------------------------------------------------------
####### 数据集 8：三个clusters low-dim, 有交叉 ########
set.seed(fixSeed)
m <- 300
dims <- 2
testData <- rbind(matrix(0.40*rnorm(m/3 * dims)+2.0,m/3,dims),
                  matrix(0.10*rnorm(m/3 * dims)+2.7,m/3,dims),
                  matrix(0.25*rnorm(m/3 * dims)+1.0,m/3,dims))
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集 9：七个clusters low-dim, 有交叉 ########
set.seed(fixSeed)
m <- 700
dims <- 2
testData <- rbind(rmvnorm(m/7, mean = rep(2.0, dims), sigma = diag(0.2^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(2.5, dims), sigma = diag(0.1^2, dims, dims)),
                  rmvnorm(m/7, mean = c(2.8, 3.0), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = c(3.0, 2.8), sigma = diag(0.05^2, dims, dims)),
                  
                  rmvnorm(m/7, mean = rep(1.5, dims), sigma = diag(0.1^2, dims, dims)),
                  rmvnorm(m/7, mean = c(1.3, 1.0), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = c(1.1, 1.4), sigma = diag(0.05^2, dims, dims))
)
plot(testData[,1], testData[,2])

#-------------------------------------------------------------------------------
####### 数据集 10：七个clusters low-dim, 交叉得更多 ########
set.seed(fixSeed)
m <- 700
dims <- 2
testData <- rbind(rmvnorm(m/7, mean = rep(2.0, dims), sigma = diag(0.2^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(2.3, dims), sigma = diag(0.1^2, dims, dims)),
                  rmvnorm(m/7, mean = c(2.5, 2.7), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = c(2.7, 2.5), sigma = diag(0.05^2, dims, dims)),
                  
                  rmvnorm(m/7, mean = rep(1.7, dims), sigma = diag(0.1^2, dims, dims)),
                  rmvnorm(m/7, mean = c(1.5, 1.3), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = c(1.3, 1.5), sigma = diag(0.05^2, dims, dims))
)
plot(testData[,1], testData[,2])


#-------------------------------------------------------------------------------
####### 数据集 10：七个clusters low-dim, 交叉得更更多，high-dim ########
set.seed(fixSeed)
m <- 700
dims <- 20
testData <- rbind(rmvnorm(m/7, mean = rep(2.0, dims), sigma = diag(0.2^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(2.3, dims), sigma = diag(0.1^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(c(2.4, 2.6), dims/2), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(c(2.6, 2.4), dims/2), sigma = diag(0.05^2, dims, dims)),
                  
                  rmvnorm(m/7, mean = rep(1.7, dims), sigma = diag(0.1^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(c(1.6, 1.4), dims/2), sigma = diag(0.05^2, dims, dims)),
                  rmvnorm(m/7, mean = rep(c(1.4, 1.6), dims/2), sigma = diag(0.05^2, dims, dims))
)
plot(testData[,1], testData[,2])

## numOfData
numOfData <- nrow(testData)

}
