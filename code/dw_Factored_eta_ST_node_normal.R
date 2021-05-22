#'@include dw_node.R
NULL

#' R6 class for Normal node.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords Normal node class
#' @field dataIds: data IDs
#' @field tssb: A TSSB_DW object
#' @field children: descentant node objects
#' @field parent: parent node object
#' @method initialize
#' @method GetChildren
Normal_DW_Factored_eta_ST <- R6::R6Class(
  classname = "Normal_DW_Factored_eta_ST",
  inherit = Node_DW_Factored_eta_ST,
  public = list(
    # Fields ------------------------------------------------------------------
    sigma = NULL,
    params = NULL,

    # Methods -----------------------------------------------------------------
    initialize = function(parent = NULL, tssb = NULL, numOfChildren = NULL, dataDims = NULL, depth = 0,  # Node_DW_Factored class
                          drift = 1,            # 这个传入的 drift 有啥用？？？
                          etaNormal = 1,
                          etaTheta  = 1,
                          priorDriftScale = 1,  
                          priorDriftShape = 1,  # drift ~ invGamma(priorDriftShape, priorDriftScale)
                          priorSigmaScale = 1,  
                          priorSigmaShape = 1,  # sigma_eps ~iid invGamma(priorSigmaShape, priorSigmavScale)
                          priorEtaNormalA = 1,
                          priorEtaNormalB = 1,  # eta_normal ~ Beta(A, B)
                          priorEtaThetaA  = 1,
                          priorEtaThetaB  = 1,  # eta_theta  ~ Beta(A, B)
                          initMean = array(0, dim=c(1, dataDims))) {
      super$initialize(parent = parent, 
                       tssb = tssb, 
                       numOfChildren = numOfChildren,
                       dataDims = dataDims,
                       depth = depth)
      # Initialize numOfChildren in Node_DW class (inherit from father class)!

      if (is.null(parent)) {
        # drift, eta_normal, eta_theta, 和对应的超参数，都只存储于根节点 root node 中！
        private$drift <- rinvgamma(dataDims, shape = priorDriftShape, scale = priorDriftScale)  # row vec
        private$etaNormal <- etaNormal
        private$etaTheta  <- etaTheta
        private$priorDriftScale <- priorDriftScale
        private$priorDriftShape <- priorDriftShape
        private$priorSigmaScale <- priorSigmaScale
        private$priorSigmaShape <- priorSigmaShape
        private$priorEtaNormalA <- priorEtaNormalA
        private$priorEtaNormalB <- priorEtaNormalB
        private$priorEtaThetaA  <- priorEtaThetaA
        private$priorEtaThetaB  <- priorEtaThetaB
        private$initMean <- initMean
        self$sigma <- rinvgamma(dataDims, shape = priorSigmaShape, scale = priorSigmaScale)  # row vec
        # self$params <- unlist(Map(function(i) {
        #                  rnorm(1, mean=initMean[i], sd = sqrt(private$drift[i]))
        #                  }, 1:dataDims))
        self$params <- rnorm(dataDims, mean = initMean, sd = sqrt(private$drift))
      } else {
        drift <- self$GetDrift()
        etaTheta  <- self$GetEtaTheta()
        self$sigma <- rinvgamma(dataDims, 
                                shape = self$GetPriorSigmaShape(), 
                                scale = self$GetPriorSigmaScale())
        # self$params <- unlist(Map(function(i) {
        #                  rnorm(1, mean=parent$params[i], sd = sqrt(etaTheta^(depth) * drift[i]))
        #                  }, 1:dataDims))
        self$params <- rnorm(dataDims, mean = parent$params, sd = sqrt(etaTheta^(self$depth) * drift))
      }
    },

    ### drift
    GetPriorDriftShape = function() {
      if (is.null(private$parent)) {
        private$priorDriftShape
      } else {
        private$parent$GetPriorDriftShape()
      }
    },

    GetPriorDriftScale = function() {
      if (is.null(private$parent)) {
        private$priorDriftScale
      } else {
        private$parent$GetPriorDriftScale()
      }
    },

    GetDrift = function() {
      if (is.null(private$parent)) {
        private$drift
      } else {
        private$parent$GetDrift()
      }
    },

    ### sigma
    GetPriorSigmaShape = function() {
      if (is.null(private$parent)) {
        private$priorSigmaShape
      } else {
        private$parent$GetPriorSigmaShape()
      }
    },

    GetPriorSigmaScale = function() {
      if (is.null(private$parent)) {
        private$priorSigmaScale
      } else {
        private$parent$GetPriorSigmaScale()
      }
    },
    
    ### eta_normal
    GetPriorEtaNormalA = function() {
      if (is.null(private$parent)) {
        private$priorEtaNormalA
      } else {
        private$parent$GetPriorEtaNormalA()
      }
    },
    
    GetPriorEtaNormalB = function() {
      if (is.null(private$parent)) {
        private$priorEtaNormalB
      } else {
        private$parent$GetPriorEtaNormalB()
      }
    },
    
    GetEtaNormal = function() {
      if (is.null(private$parent)) {
        private$etaNormal
      } else {
        private$parent$GetEtaNormal()
      }
    },
    
    ### eta_theta
    GetPriorEtaThetaA = function() {
      if (is.null(private$parent)) {
        private$priorEtaThetaA
      } else {
        private$parent$GetPriorEtaThetaA()
      }
    },
    
    GetPriorEtaThetaB = function() {
      if (is.null(private$parent)) {
        private$priorEtaThetaB
      } else {
        private$parent$GetPriorEtaThetaB()
      }
    },
    
    GetEtaTheta = function() {
      if (is.null(private$parent)) {
        private$etaTheta
      } else {
        private$parent$GetEtaTheta()
      }
    },

    # Log probability
    GetLogProb = function(x) {
      # 没有datum时，传入的 x 是一个空矩阵(matrix array类型)，dnorm() 返回 numeric(0)，取sum返回0
      etaNormal <- self$GetEtaNormal()
      depth <- self$depth
      
      if (!is.null(dim(x))) {
        sum(apply(x, MARGIN = 1, function(datum) {dnorm(datum, mean = self$params, sd = sqrt(etaNormal^(depth) * self$sigma), log = T)} ))
      } else {
        sum(dnorm(x, mean = self$params, sd=sqrt(etaNormal^(depth) * self$sigma), log = T))
      }
    },

    GetNodeLogProb = function() {
      self$GetLogProb(self$GetData())
    },
    
    
    # param_log_prior
    GetParamLogPrior = function() {
      etaTheta <- self$GetEtaTheta()
      mean <- if (is.null(private$parent)) private$initMean else private$parent$params
      dnorm(self$params, mean = mean, sd = sqrt(etaTheta^(self$depth) * self$GetDrift()))
    },
    
    # # Non-log probability
    # GetProb = function(x) {
    #   sum(dmvnorm(x, self$params, self$sigma, log = FALSE))
    # },
    # 
    # GetNodeProb = function() {
    #   self$GetProb(self$GetData())
    # },

    ### Update params (theta and sigma)
    # Input:
    #    keep logical vector of self$children whether needs to be updated
    #   （keep 也可以不传入，通过 HasData 函数达到相同的效果）
    #   （但如果不传入 keep，而是用 HasData 函数，就导致每个 node 都需要从它的所有children
    #     递归搜索后返回 T/F，如果节点很多，也会很慢。而 keep 是在 KeepTree 函数作用后一次性生成好的)
    ResampleParams = function(keep) {
      nodeData <- self$GetData()
      drift <- self$GetDrift()
      etaNormal <- self$GetEtaNormal()
      etaTheta  <- self$GetEtaTheta()
      depth <- self$depth
      
      numOfData <- self$GetNumOfLocalData()  # Number of local data at this node (self)
      priorSigmaScale <- self$GetPriorSigmaScale()
      priorSigmaShape <- self$GetPriorSigmaShape()
      numOfChildren <- sum(keep)
      # actChildrenInd <- Filter(function(i) private$children[i]$HasData(), 1:self$numOfChildren)
      # numOfActChildren <- length(actChildrenInd)

      # Data
      if (numOfData == 0) {
        dataMean = 0
      } else if (numOfData == 1){
        dataMean = nodeData
      } else {
        dataMean = colMeans(nodeData)
      }

      # Parent param
      if (is.null(private$parent)) {
        parentParams <- private$initMean
      } else {
        parentParams <- private$parent$params
      }

      # Children params
      # numOfChildren = W
      if (numOfChildren ==0 ) {  # numOfActChildren
        childParamsMean = 0
      } else {
        childParams <- Reduce(
          rbind,
          Map(function(x) {x$params}, private$children[keep]),
          # Map(function(x) {x$params}, private$children[actChildrenInd]),
          c()
        )
        childParamsMean <- colMeans(childParams)
      }

      # Construct prior for node mean
      priorParamsMean <- (etaTheta*parentParams + numOfChildren*childParamsMean)/(numOfChildren + etaTheta)
      priorParamsCov <- etaTheta^(depth+1)*drift / (numOfChildren + etaTheta)  # drift 每个分量都除以

      # Posterior for node mean
      if (numOfData == 0) {
        postParamsMean <- priorParamsMean
        postParamsCov <- priorParamsCov
      } else {
        postParamsCov <- (priorParamsCov^(-1) + numOfData/(etaNormal^(depth) * self$sigma))^(-1)  # 此处所有运算都是按元素求的
        postParamsMean <- (priorParamsMean/priorParamsCov +
                             numOfData*dataMean/(etaNormal^(depth) * self$sigma)) * postParamsCov
      }
      self$params <- rnorm(n = self$dataDims,
                           mean = postParamsMean,
                           sd = sqrt(postParamsCov))

      # Posterior for node covarirance
      # 注意：
      #   riwish 中的 S 如果不是对称矩阵（同时 chol 函数不报错）
      #   那么抽取出的 sigma 也不对称！！！
      if (numOfData > 0) {  
        # 如果 numOfData = 0，不需要从 prior 中更新吗？？？？
        postSigmaShape <- priorSigmaShape + numOfData/2
        
        dataSqu <- if (numOfData == 1) {
          (nodeData - self$params)^2
        } else {
          rowSums(apply(nodeData, MARGIN = 1, function(x) {(x - self$params)^2} ))
        }
        postSigmaScale <- priorSigmaScale + dataSqu / (2*etaNormal^(depth))
        self$sigma <- rinvgamma(n = self$dataDims, shape = postSigmaShape, scale = postSigmaScale) 
      }
      invisible(self)
    },

    ### Update drift
    # Finite nodes: No need to use SliceSampler !!
    # maxDepth (in tssb) can get by: self$tssb$maxDepth
    ResampleHyperDrift = function(verbose = F) {
      if (verbose) cat("Resample drift Start...\n")
      priorDriftScale <- self$GetPriorDriftScale()
      priorDriftShape <- self$GetPriorDriftShape()
      etaTheta <- self$GetEtaTheta()
      # depth <- self$depth

      if(!is.null(private$parent)) {
        stop("Can only update hypers from root!")  # 因为只有root node，才能调用 private$drift！
      }
      
      Descend <- function(root, depth = 1) {
        llh <- 0
        numOfActNodes <- 0
        actChildren <- Filter(function(x) {x$HasData()}, root$GetChildren())

        for (i in seq_along(actChildren)) {
          # 如果用 1:numOfActChildren，则会出现 i=1，i=0！
          numOfActNodes <- numOfActNodes + 1
          child <- actChildren[[i]]
          res <- Descend(child, depth = depth+1)
          llh <- llh + (child$params - root$params)^2 / etaTheta^(depth)  # 循环中先不除2！最后除一次即可
          llh <- llh + res$llh
          numOfActNodes <- numOfActNodes + res$num
        }
        return(list(llh = llh, num = numOfActNodes))
      }
      
      res <- Descend(self)
      postDriftShape <- priorDriftShape + (res$num + 1) / 2  # +1 for root node!
      postDriftScale <- priorDriftScale + (res$llh + (self$params - private$initMean)^2) / 2
      # root 和 initMean 之间也有 drift，但是没有 etaTheta！
      
      private$drift <- rinvgamma(n = self$dataDims, shape = postDriftShape, scale = postDriftScale)
      if (verbose) cat("Resample drift Finish!\n")
    },
    
    
    # Update two Eta
    ResampleHyperEta = function(sampleEtaNormal = T,
                                sampleEtaTheta  = T,
                                verbose = F) {
      if (verbose) cat("Resample two eta Start...\n")
      if(!is.null(private$parent)) {
        stop("Can only update hypers from root!")
      }
      priorEtaNormalA <- self$GetPriorEtaNormalA()
      priorEtaNormalB <- self$GetPriorEtaNormalB()
      priorEtaThetaA  <- self$GetPriorEtaThetaA()
      priorEtaThetaB  <- self$GetPriorEtaThetaB()
      drift <- self$GetDrift()
      
      ComputeEtaNormalLlh <- function(etaNormal) {
        if (etaNormal <= 0 || etaNormal >= 1) {
          return(-Inf)
        }
        res <- self$tssb$GetLogCompleteDataLikelihood(flagRoot=F)
        return(res$ll + dbeta(etaNormal, shape1 = priorEtaNormalA, shape2 = priorEtaNormalB, log = TRUE))
      }
      
      ComputeEtaThetaLlh <- function(etaTheta) {
        if (etaTheta <= 0 || etaTheta >= 1) {
          return(-Inf)
        }
        Descend <- function(root, depth = 1) {
          llh <- 0
          actChildren <- Filter(function(x) {x$HasData()}, root$GetChildren())
          for (i in seq_along(actChildren)) {
            child <- actChildren[[i]]
            res <- Descend(child, depth=depth+1)
            llh <- llh + sum(dnorm(child$params, mean = root$params, sd = sqrt(etaTheta^(depth) * drift), log = TRUE))
            llh <- llh + res
          }
          return(llh)
        }
        return(Descend(self) + dbeta(etaTheta, shape1 = priorEtaThetaA, shape2 = priorEtaThetaB, log = TRUE))
      }
      
      if (sampleEtaNormal) {
        private$etaNormal <- SliceSampler(private$etaNormal, ComputeEtaNormalLlh, verbose = verbose)
        if (verbose) cat("Resample eta_normal Finish!\n")
      }
      if (sampleEtaTheta) {
        private$etaTheta  <- SliceSampler(private$etaTheta,  ComputeEtaThetaLlh,  verbose = verbose)
        if (verbose) cat("Resample eta_theta Finish!\n")
      }

    }
    
  ),
  private = list(
    drift = NA,
    etaNormal = NA,
    etaTheta = NA,
    priorDriftScale = NA,
    priorDriftShape = NA,
    priorSigmaScale = NA,
    priorSigmaShape = NA,
    priorEtaNormalA = NA,
    priorEtaNormalB = NA,
    priorEtaThetaA  = NA,
    priorEtaThetaB  = NA,
    initMean = NA
  )
)