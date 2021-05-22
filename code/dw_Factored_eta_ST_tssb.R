#'TSSB_DW is a R6 object of tree-structured stick breaking
#' process
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @field data
#' @field dpAlpha
#' @field dpGamma
#' @field dpLambda
#' @field minDepth
#' @field maxDepth
#' @field root
#' @field assignments
#' @field numOfData
#' @field minDpAlpha
#' @field maxDpAlpha
#' @field minDpLambda
#' @field minDpLambda
#' @field minDpGamma
#' @field maxDpGamma
#' @method new TSSB_DW constuctor
#' @method FindNode Find node in tree, generate new node as needed
#' @method GetMixture Compute the mixing weights and get the corresponding nodes
#' @method ConvertTssbToIgraph Convert the tssb root list to an igraph object
#' @method CullTree Remove empty leaf nodes
#' @method GetLogMarginalDataLikelihood Compute the log marginal data likelihood
TSSB_DW_Factored_eta_ST <- R6::R6Class(
  classname = "TSSB_DW_Factored_eta_ST",
  public = list(
    data = NULL,   # must be matrix array type! (if one dimension) !!!!
    # Felids ------------------------------------------------------------------
    dpAlpha = NA,
    dpGamma = NA,
    dpLambda = NA,
    minDepth = NA,
    maxDepth = NA,
    root = NULL,
    assignments = NULL,  # List of nodes
    numOfData = NA,
    minDpAlpha = NA,
    maxDpAlpha = NA,
    minDpLambda = NA,
    maxDpLambda = NA,
    minDpGamma = NA,
    maxDpGamma = NA,
    
    # new added!!!
    maxWidth = NA, 

    # Initialize --------------------------------------------------------------
    initialize = function(rootNode = emptyenv(),
                          dpAlpha = 1,
                          dpGamma = 1,
                          dpLambda = 1,
                          data = NULL,
                          minDepth = 0,
                          maxDepth = 5,
                          maxWidth = 3,  # i.e. W new added!!!
                          minDpAlpha = 1e-6, maxDpAlpha = 50,
                          minDpLambda = 1e-6, maxDpLambda = 1,
                          minDpGamma = 1e-6, maxDpGamma = 50,
                          Flag.onlyTree = T
                          ) {
      if (!is.environment(rootNode) ||
            identical(rootNode, emptyenv()) ||
            !"Node_DW_Factored_eta_ST" %in% class(rootNode)) {
        stop("Root node must be specified")
      }

      self$data    <- data  # data 只记录在 tssb 中，node 记录的是 dataIds，注意区分！
      self$dpAlpha <- dpAlpha
      self$dpGamma <- dpGamma
      self$minDepth <- minDepth
      self$maxDepth <- maxDepth
      self$maxWidth <- maxWidth  # new added!!!
      self$dpLambda <- dpLambda
      self$minDpAlpha <- minDpAlpha
      self$maxDpAlpha <- maxDepth
      self$minDpLambda <- minDpLambda
      self$maxDpLambda <- maxDpLambda
      self$minDpGamma <- minDpGamma
      self$maxDpGamma <- maxDpGamma
      self$numOfData <- if (is.null(data)) 0 else nrow(self$data)
      self$root <- list(
        node     = rootNode,
        path     = c(),    # 其实 path 在整个过程中并没有用到，只是查阅方便了
        keep     = c(),
        main     = if (self$minDepth == 0) BoundBeta(1, 1, self$dpAlpha) else 0,
        sticks   = c(),
        children = list())
      rootNode$tssb <- self
      rootNode$numOfChildren <- maxWidth  # new added !!!

      # Initialize the whole tree and data
      cat("Initialization: D =", self$maxDepth,"and W =", self$maxWidth,"\n")
      self$InitWholeTree()
      if (!is.null(data)) {
        self$InitWholeData()
        self$KeepTree(Flag.onlyTree)
      }
      cat("Initialization Is Over!\n")
    }, # end of initialize
      
    # Other Methods -----------------------------------------------------------
    
    # Initialize the whole tree
    InitWholeTree = function(){
      Descend <- function(root, depth = 0) {
        if (depth == self$maxDepth) return(root = root)
        for (i in 1:(self$maxWidth)) {
          newStick <- if (i != self$maxWidth) {BoundBeta(1, 1, self$dpGamma)} else {1.0} 
          ### 此处多次 if 判断，可能耗时？？？
          root$sticks <- c(root$sticks, newStick)
          root$children <- c(root$children,
                             list(
                               list(
                                 node = root$node$Spawn(),
                                 path = c(root$path, i),
                                 keep = c(),
                                 main = if (self$minDepth > (depth+1)) {0.0
                                 } else if (self$maxDepth > (depth+1)) {
                                   BoundBeta(1, 1, (self$dpLambda^(depth+1))*self$dpAlpha)
                                 } else {1.0},
                                  sticks = NULL,
                                 children = NULL
                               )
                             )
                           )
        }
        # root$children <- Map(Descend, root$children, depth = depth + 1)
        # 上面的代码可能速度更快？？？
        for (i in 1:(self$maxWidth)) {
          if (depth == 0) cat("No.", i, "child of root is initializing...\n")
          root$children[[i]] <- Descend(root$children[[i]], depth = depth+1)
        }
        # 将更新后的 children 重新赋值给 root$children
        return(root = root)
      }
      self$root <- Descend(self$root)
    },

    # Initialize the data in the tree
    # C_i | Pi ~ MN(Pi)
    # 由于在 truncation 版本中，datum 属于哪个 node 是从 multinoulli 分布直接抽到的，
    # 并不是 Top-down 的搜索，因此无法像 FindNode() 函数中那样，找到之后再 Down-top 的方向更新 root，
    # 而是在 抽到 C_i 之后，先 Top-down 找到这个点，然后数据赋值，再 Down-top 更新 root。
    # InitWholeData = function(){
    #   res <- self$GetMixture()
    #   # # return:
    #   # # rootList, weight
    #   # Indices <- rmultinom(self$numOfData, 1, res$weight)
    #   weight <- res$weight
    #   numOfData <- sum(weight != 0)
    #   weight[weight!=0] <- 1/numOfData
    #   
    #   
    #   ##### 初始化时均匀分到所有“有权重”的nodes上
    #   Indices <- rmultinom(self$numOfData, 1, weight)   
    #   # len(weight) * numOfData array matrix
    #   # Each column is a multinom sample
    #   rootSeq <- apply(Indices, MARGIN = 2, function(x) {res$rootList[[which(x == 1)]] } )
    #   # return a list of N roots (lists)
    #   self$assignments <- Map(function(x) {x$node}, rootSeq)
    # 
    #   # Add data and Update roots
    #   # Add data (This procedure DO NOT need to update root!)
    #   # (Only when produce new children should we update root!)
    #   Map(function(n) {self$assignments[[n]]$AddDatum(n)}, 1:self$numOfData)
    #   # Descend <- function(root, path, n, depth = 0) {
    #   # # path 和 n 是否传入，哪个更快？？？
    #   #   if (depth == length(path)) {
    #   #     root$node$AddDatum(n)
    #   #     return(root = root)
    #   #   }
    #   #   root$children[[ path[depth+1] ]] <- Descend(root$children[[ path[depth+1] ]], path, n, depth=depth+1)
    #   #   return(root = root)
    #   # }
    #   # 
    #   # for (n in 1:self$numOfData) {
    #   #   path <- rootSeq[[n]]$path
    #   #   self$root <- Descend(self$root, path, n)
    #   # }
    # },
    
    # Initialize all data in the nodes evenly !!
    InitWholeData = function(){
      res <- self$GetMixture()
      rootList <- res$rootList
      weight <- res$weight
      numOfData  <- self$numOfData
      dataInd <- sample(rep_len(which(weight != 0), length.out = numOfData))  # random permutation

      self$assignments <- Map(function(i) {rootList[[i]]$node}, dataInd)
      Map(function(n) {self$assignments[[n]]$AddDatum(n)}, 1:numOfData)
    },
      
    # GetMixture function
    # Return rootList (list of roots (lists)), not the nodes vector!!!
    # 但作用和 Yuan 代码中的 nodes 的作用相同
    ###########
    # 注意：
    #    1. 使用 "log-sum-exp" trick ! ! !
    #    2. NOT input datum
    ###########
    GetMixture = function() {
      Descend <- function(root, logMass) {
        rootList <- list(root)
        logWeight <- logMass + log(root$main)
        # if (!is.null(datum)) {
        #   logWeight <- logWeight + root$node$GetLogProb(datum)
        # } 
        edges <- SticksToEdges(root$sticks)
        weights <- diff(c(0, edges))

        for (i in seq_along(root$children)) {
        # 改成 1:self$numOfChildren 加 depth 判断可能更快？？？
          child <- root$children[[i]]
          res <- Descend(child, logMass + log(1.0-root$main) + log(weights[i]))
          rootList <- c(rootList, res$rootList)
          logWeight <- c(logWeight, res$logWeight)
        }
        return(list(rootList = rootList, logWeight = logWeight))
      }
      
      res <- Descend(self$root, 0.0)  # 此时，res$logWeight 存储的是 log(权重)
      maxLogWeight <- max(res$logWeight)
      weight <- exp(res$logWeight - maxLogWeight) / sum(exp(res$logWeight - maxLogWeight))
      return(list(rootList = res$rootList, weight = weight))
    },
    
    
    # Compute log-likelihood of datum X_i at every node: 
    #   p(X_i | theta_eps, sigma^2_eps), for ALL eps (including inactive nodes)
    GetDatumCurrentLogLikelihood = function(datum) {
      Descend <- function(root) {
        rootList <- list(root)
        llhList  <- root$node$GetLogProb(datum)
        
        for (i in seq_along(root$children)) {
          # 改成 1:self$numOfChildren 加 depth 判断可能更快？？？
          child <- root$children[[i]]
          res <- Descend(child)
          rootList <- c(rootList, res$rootList)
          llhList  <- c(llhList,  res$llhList)
        }
        return(list(llhList = llhList, rootList = rootList))
      }
      return(Descend(self$root))
    },
    

    ConvertTssbToIgraph = function() {
      edges <- SticksToEdges(self$root$sticks)
      weights <- diff(c(0, edges))
      g <- igraph::graph.empty(directed = TRUE)
      g <- g + igraph::vertex(name = "X",
                              size = self$root$node$GetNumOfLocalData())

      Descend <- function(root, name, mass, g) {
        if (length(root$sticks) < 1){
          return(list(total = mass, g=g))
        } else {
          total <- 0
          edges <- SticksToEdges(root$sticks)
          weights <- diff(c(0, edges))

          for (i in 1:length(root$children)) {
            child <- root$children[[i]]
            childName <- paste(name, i, sep = "-")
            childMass <- mass * weights[i] * child$main
            g <- g + igraph::vertex(name = childName,
                                    size = child$node$GetNumOfLocalData())
            g <- g + igraph::edge(name, childName,
                                  Value = child$node$GetNumOfLocalData())
            tmp <- Descend(child,
                           childName,
                           mass*weights[i]*(1.0 - child$main),
                           g)
            g <- tmp$g
            total = total + childMass + tmp$total
          }
          return(list(total=total, g=g))
        }
      }
    res = Descend(self$root, "X", 1, g)
    return(res)
    },

    # Method: KeepTree
    #   Find the nodes that in the tree structure (corresponding keep = T/F)
    #   Flag.onlyTree = T: update only the active nodes
    #                   F: update the whole DW tree
    KeepTree = function(Flag.onlyTree = T) {
      Descend <- function(root) {
        res <- unlist(Map(Descend, root$children), recursive = F, use.names = F)
        # 返回的是 2 * length(children) 个元素的list，其中
        #    奇数项为 counts，偶数项为 root

        if (length(res) == 0) {  # leaf node
          return(list(counts = root$node$GetNumOfLocalData(), root = root))
        }
        counts <- unlist(res[seq(1, length(res), by = 2)])  # 1, 3, 5
        root$children <- res[seq(2, length(res), by = 2)]   # 2, 4, 6  # 更新 root
        # keep <- which(counts != 0)  # 记录含有data 的nodes
        
        root$keep <- if (Flag.onlyTree) (counts != 0) else rep(TRUE, self$maxWidth)

        return(list(counts = sum(counts) + root$node$GetNumOfLocalData(),
                    root = root))  
        # 将 counts 累积起来返回，这样即使中间节点没有数据，若其child有数据，也可以保证不被kill
      }
      res <- Descend(self$root)
      self$root <- res$root
      invisible(self)
    },
    
    # Method: CullTree
    #   cull the tree to kill the nodes not in the tree structure
    CullTree = function() {
      Descend <- function(root) {
        res <- unlist(Map(Descend, root$children), recursive = F, use.names = F)
        
        if (length(res) == 0) {
          return(list(counts = root$node$GetNumOfLocalData(), root = root))
        }
        counts <- unlist(res[seq(1, length(res), by = 2)])
        root$children <- res[seq(2, length(res), by = 2)]
        keep <- which(counts != 0)
        
        sapply(root$children[which(counts == 0)], function(x) x$node$Kill())
        
        if (length(keep) == 0) {
          root$sticks <- NULL
          root$children <- NULL
        } else {
          root$sticks <- root$sticks[keep]
          root$children <- root$children[keep]
        }
        return(list(counts = sum(counts) + root$node$GetNumOfLocalData(),
                    root = root))
      }
      res <- Descend(self$root)
      self$root <- res$root
      invisible(self)
    },

    GetLogMarginalDataLikelihood = function() {
      res <- self$GetMixture()
      weights <- res$weight
      roots <- res$rootList

      ll <- Reduce(
        sum,
        Map(function(i) {
          node <- roots[[i]]$node
          if (node$GetNumOfLocalData()) {
            node$GetNumOfLocalData()*weights[i] + node$GetNodeLogProb()
            # 为什么要加第一项？？？
            # node$GetNodeLogProb()
          } else {
            0
          }
        },
        seq_along(weights)
        ),
        0
      )
      return(list(ll=ll, ww = weights, nn = roots))
    },
    
    GetLogCompleteDataLikelihood = function(flagRoot=T) {
      res <- self$GetMixture()
      weights <- res$weight
      roots <- res$rootList
      seqWeights <- if (flagRoot==T) seq_along(weights) else seq_along(weights)[-1]
      
      ll <- Reduce(
        sum,
        Map(function(i) {
          node <- roots[[i]]$node
          if (node$GetNumOfLocalData()) {node$GetNodeLogProb()} else {0.0}
        },
        seqWeights
        ),
        0
      )
      return(list(ll=ll, ww = weights, nn = roots))
    },
    
    # Compute the unnormalized posterior
    # likelihood * prior
    GetUnnormalizedPost = function() {
      res <- self$GetMixture()
      weights <- res$weight
      roots <- res$rootList
      
      ll <- Reduce(
        sum,
        Map(function(i) {
          node <- roots[[i]]$node
          if (node$GetNumOfLocalData()) {
            node$GetNumOfLocalData()*log(weights[i]) + 
              node$GetNodeLogProb() +
              node$GetParamLogPrior()
          } else {
            0
          }
        },
        seq_along(weights)
        ),
        0
      )
      return(list(ll=ll, ww = weights, nn = roots))
    }
  )
)



