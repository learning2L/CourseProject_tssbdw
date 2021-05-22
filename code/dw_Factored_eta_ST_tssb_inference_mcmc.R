#'@include tssb.R
NULL

#' R6 class for inference via MCMC for TSSB.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @method ResampleAssignments Sample data assignments with Adams's slice-retrospective sampler
#' @method ResampleSticks Sample sticks
#' @method ResampleStickOrders Sample stick orders
#' @method ResampleNodeParameters Sample node parameters
#' @method ResampleHypers Sample TSSB hyper parameters

TssbMCMC_DW_Factored_eta_ST <- R6::R6Class(
  classname = "TssbMCMC_DW_Factored_eta_ST",
  inherit = TSSB_DW_Factored_eta_ST,
  public = list(
    initialize = function(...) {
      super$initialize(...)
    },  # end of initialize
    
    
    
    
    # Compute the post weights
    # using log-sum-exp trick !
    ComputePostWeight = function(weight, llhList) {
      logPostWeight <- log(weight) + llhList
      maxElem <- max(logPostWeight)
      postWeight <-exp(logPostWeight - maxElem) / sum(exp(logPostWeight - maxElem))
      return(postWeight)
    },
    
    
    
    # 关键：如何将原数据位置 remove 掉！
    # 解答：
    #      1. 直接用 self$assignment 来 remove data 即可，不需要 update root！
    #      2. 这是因为 R6 class 直接赋值，只是传递了引用，并没有 clone！
    #      3. 但如果改变了 root 中的其他变量值(即改变了tssb)，就需要更新 root！（如更新sticks）
    # 注意：
    #      当节点很多，或者likelihood接近于0时，GetMixture(datum) 要使用 "log-sum-exp" 技巧！
    #      否则 weights 可能会出现 0！
    ResampleAssignments = function(verbose=FALSE) {
      if (verbose) cat("Resample Assignment Start...\n")
      res <- self$GetMixture()
      priorWeight <- res$weight
      rootList <- res$rootList
      
      for (n in 1:self$numOfData) {
        # 此 for 循环可以考虑用 parallel computing！
        resDatum <- self$GetDatumCurrentLogLikelihood(self$data[n,])
        weight <- self$ComputePostWeight(priorWeight, resDatum$llhList)
        Index <- rmultinom(1, 1, weight)
        NewRoot <- unlist(rootList[as.logical(Index)], recursive = F)
        NewAssignment <- NewRoot$node
        
        if (!identical(NewAssignment, self$assignments[[n]])) {
          newNode <- NewAssignment
          self$assignments[[n]]$RemoveDatum(n)
          newNode$AddDatum(n)
          self$assignments[[n]] <- newNode
        }
      }
      if (verbose) cat("Resample Assignment Finish!\n")
      invisible(self)
    },

    ResampleSticks = function(verbose=F) {
      if (verbose) cat("Resample Sticks Start...\n")
      Descend <- function(root, depth=0) {
        dataDown <- 0
        # indices <- seq_along(root$children)
        # deIndices <- sort(indices, decreasing = T)
        actInd <- Filter(function(n) {root$keep[n]}, 1:self$maxWidth)
        actInd.dec <- sort(actInd, decreasing = T)

        # for (i in deIndices) {  # 此处 decreasing=T 很重要！
        for (i in actInd.dec) { 
          res <- Descend(root$children[[i]], depth+1)
          childData <- res$nodeData
          root$children[[i]] <- res$root
          postAlpha <- 1 + childData
          postBeta <- self$dpGamma + dataDown
          root$sticks[i] <- if (i < self$maxWidth) BoundBeta(1, postAlpha, postBeta) else {1.0}
          # psi_{eW} 恒为 1
          dataDown <- dataDown + childData
        }

        dataHere <- root$node$GetNumOfLocalData()
        postAlpha <- 1 + dataHere
        postBeta <- (self$dpLambda^depth) * self$dpAlpha + dataDown
        root$main <- if (self$minDepth > depth) {0.0
        } else if (self$maxDepth > depth) {
          BoundBeta(1, postAlpha, postBeta)
        } else {1.0}
        
        # nu_eps 恒为 0，当 |eps| < minDepth
        # nu_eps 恒为 1，当 |eps| = maxDepth
        return(list(nodeData = dataHere + dataDown, root = root))
      }

      self$root = Descend(self$root)$root
      if (verbose) cat("Resample Sticks Finish!\n")
      invisible(self)
    },

    
    # Down-Top update!
    # Only update the nodes in the tree structure (keep = T, named "active nodes")
    ResampleNodeParameters = function(verbose=FALSE) {
      if (verbose) cat("Resample Node Params Start...\n")
      Descend <- function(root) {
        # 无法从外部访问 Node class 的私有成员 children
        # 所以无法调用 root$node$children
        # 但可以用 root$node$GetChildren 的方式取出children，所以可能不需要 keep变量
        Map(Descend, root$children[root$keep])
        root$node$ResampleParams(root$keep)
        # 只更新 keep = T 的nodes！！
        # 第二行其实可以不用传入 root$keep，直接通过调用 node 类的 HasData 函数也可以达到相同目的
        
        return(root=root)
      }
      self$root <- Descend(self$root)
      if (verbose) cat("Resample Node Params Finish!\n")
      invisible(self)
    },

    ResampleHypers = function(sampleDpAlpha = T,
                              sampleDpLambda = T,
                              sampleDpGamma = T,
                              verbose = F) {
    # 用 Slice sampler 更新这三个超参数时要注意：
    # 因为在 TSSB-DW 中，最后一个 child 的 psi 恒为1，最后一层 children 的 main 恒为1
    # 必须要考虑这一点，否则 llh 总是 -Inf
      if (verbose) cat("Resample Hypers Start...\n")

      # Compute log-likelihood of v_e
      # \prod_{e} Be(v_e | 1, labmda^{depth} * alpha)
      # Notice: v_e = 1, for all |e| = maxDepth. Should not compute likelihood!
      ComputeDpMainLlh <- function(params,
                                   sample = c("alpha", "lambda"),
                                   fixParams = NULL) {
        if (sample == "alpha") {
          dpAlpha <- params
          dpLambda <- fixParams
        } else if (sample == "lambda") {
          dpAlpha <- fixParams
          dpLambda <- params
        } else {
          stop("sample has to be either 'alpha' or 'lambda' ")
        }

        # minAlpha < alpha < maxAlpha
        if (dpAlpha < self$minDpAlpha || dpAlpha > self$maxDpAlpha) {
          return(-Inf)
        }
        # minLambda < lambda < maxLambda
        if (dpLambda < self$minDpLambda || dpLambda > self$maxDpLambda) {
          return(-Inf)
        }
        Descend <- function(root, depth = 0) {
          llh <- if (self$minDepth <= depth) {
            dbeta(root$main, 1, dpLambda^depth*dpAlpha, log = TRUE)
          } else {0}
          if (depth == (self$maxDepth-1)) return(list(llh = llh, root = root)) # 避免出现 which(NULL)

          for (i in which(root$keep)) {  
            # only active nodes !
            res <- Descend(root$children[[i]], depth + 1)
            llh <- llh + res$llh
            root$children[[i]] <- res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }

      # Compute log-likelihood of psi_e
      # \prod_{e} Be(psi_e | 1, gamma)
      # Notice: psi_{ej} = 1, for all j = maxWidth. Should not compute likelihood!
      ComputeDpBranchLlh <- function(dpGamma) {
        if (dpGamma < self$minDpGamma || dpGamma > self$maxDpGamma) {
          return(-Inf)
        }
        Descend <- function(root) {
          llh <- 0
          # indices <- seq_along(root$children)
          
          actInd <- Filter(function(n) {root$keep[n]}, 1:self$maxWidth)
          # 若 keep 全是 F，则 actInd 返回 integer(0)，下面 for循环不执行，即停止再向下递归
          # 若 keep 为 NULL，则 actInd 也返回 integer(0)
          # 避免了 which(NULL) 报错的情况！（which必须传入逻辑值）
          
          for (i in actInd) {
            if (i < self$maxWidth) {
              llh <- llh + dbeta(root$sticks[i], 1, dpGamma, log = TRUE)
            }
            res <- Descend(root$children[[i]])
            llh <- llh + res$llh
            root$children[[i]] <- res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }

      if (sampleDpAlpha) {
        self$dpAlpha <- SliceSampler(self$dpAlpha,
                                     ComputeDpMainLlh,
                                     sample = "alpha",
                                     fixParams = self$dpLambda,
                                     verbose = verbose)
        if (verbose) cat("  Resample alpha  Finish!\n")
      }
      if (sampleDpLambda) {
        self$dpLambda <- SliceSampler(self$dpLambda,
                                      ComputeDpMainLlh,
                                      sample = "lambda",
                                      fixParams = self$dpAlpha,
                                      verbose = verbose)
        if (verbose) cat("  Resample lambda Finish!\n")

      }
      if (sampleDpGamma) {
        self$dpGamma <- SliceSampler(self$dpGamma, ComputeDpBranchLlh,
                                     verbose = verbose)
        if (verbose) cat("  Resample gamma  Finish!\n")
      }
      invisible(self)
    },
    
    
    # Search Tree
    # [Yuan2015] swap-move
    SearchTree = function(swapMain = T, verbose = F) {
      res <- self$GetUnnormalizedPost()
      post <- rexp(1) + res$ll
      # 此处 rexp 等价于 log(unif)，为什么有这么一项？？？？
      # 是基于 slice sampler 的吗？？？所以是 MH + slice sampler 框架？？
      weight <- res$ww
      rootList <- res$nn
      len <- length(rootList)
      
      rootHasDataBool <- unlist(Map(function(i) {length(rootList[[i]]$node$dataIds) > 0}, 1:len))  # 记录的是有数据的节点，而不是active nodes！
      rootHasData <- which(rootHasDataBool)  # 包括没有数据但是在 tree-struct 中的nodes，和不在 tree-struct 的nodes
      rootHasNotData <- which(!rootHasDataBool)
      
      if (len > 1) {
        if (length(rootHasNotData) > 0) {
          nodeAnum <- sample(rootHasData, 1)
          nodeBnum <- sample(rootHasNotData, 1)
        } else {  # ALL nodes have data!
          nodeAnum <- sample(1:len, 1)
          nodeBnum <- sample(1:len, 1)
          while (nodeAnum == nodeBnum) nodeBnum <- sample(1:len, 1)
        }
        
        SwapNodes <- function(nodeAnum, nodeBnum) {
          
          
          rootA <- rootList[[nodeAnum]]
          rootB <- rootList[[nodeBnum]]
          
          paramsA <- rootA$node$params
          sigmaA  <- rootA$node$sigma
          dataIdsA <- rootA$node$dataIds 
          
          rootA$node$params <- rootB$node$params
          rootA$node$sigma  <- rootB$node$sigma
          
          Map(function(i) {rootA$node$RemoveDatum(i)}, rootA$node$dataIds)
          Map(function(i) {rootA$node$AddDatum(i)
                           self$assignments[[i]] <- rootA$node}, rootB$node$dataIds)
          
          rootB$node$params <- paramsA
          rootB$node$sigma  <- sigmaA
          
          Map(function(i) {rootB$node$RemoveDatum(i)}, rootB$node$dataIds)
          Map(function(i) {rootB$node$AddDatum(i)
                           self$assignments[[i]] <- rootB$node}, dataIdsA)
          
          # swap mains
          if (swapMain) {
            Descend <- function(root, nodeNum) {
              if (nodeNum == nodeAnum) root$main <- mainB
              if (nodeNum == nodeBnum) root$main <- mainA
              children <- root$children
              for (i in seq_along(children)) {
                nodeNum <- nodeNum + 1
                res <- Descend(children[[i]], nodeNum)
                root$children[[i]] <- res$root
                nodeNum <- res$nodeNum
              }
              return(list(root=root, nodeNum=nodeNum))
            }
            
            if (rootA$main != 1 && rootB$main != 1) {
              mainA <- rootA$main
              mainB <- rootB$main
              self$root <- Descend(root=self$root, nodeNum=1)$root
            }
          }
        }
        
        SwapNodes(nodeAnum, nodeBnum)
        self$ResampleSticks()
        post_new <- self$GetUnnormalizedPost()$ll
        if (post_new < post) {
          SwapNodes(nodeAnum, nodeBnum)
          self$ResampleSticks()
        } else {if (verbose) cat("Successful swap!\n")}
      }
    },
    
    
    # Search Tree
    # Yuan2015 思路完全复刻
    SearchTreeYuan = function(swapMain = T, verbose = F) {
      res <- self$GetUnnormalizedPost()
      post <- rexp(1) + res$ll
      # 此处 rexp 等价于 log(unif)，为什么有这么一项？？？？
      # 是基于 slice sampler 的吗？？？所以是 MH + slice sampler 框架？？
      weight <- res$ww
      rootList <- res$nn
      
      actRoots <- Filter(function(x) {x$node$HasData()}, rootList)
      len <- length(actRoots)
      nodeAnum <- sample(1:len, 1)
      nodeBnum <- sample(1:len, 1)
      while (nodeAnum == nodeBnum) nodeBnum <- sample(1:len, 1)
      
      if (len > 1) {
        SwapNodes <- function(nodeAnum, nodeBnum) {
          # swap nodeParams and dataIds
          rootA <- actRoots[[nodeAnum]]
          rootB <- actRoots[[nodeBnum]]
          
          paramsA <- rootA$node$params
          sigmaA  <- rootA$node$sigma
          dataIdsA <- rootA$node$dataIds 
          
          rootA$node$params <- rootB$node$params
          rootA$node$sigma  <- rootB$node$sigma
          
          Map(function(i) {rootA$node$RemoveDatum(i)}, rootA$node$dataIds)
          Map(function(i) {rootA$node$AddDatum(i)
            self$assignments[[i]] <- rootA$node}, rootB$node$dataIds)
          
          rootB$node$params <- paramsA
          rootB$node$sigma  <- sigmaA
          
          Map(function(i) {rootB$node$RemoveDatum(i)}, rootB$node$dataIds)
          Map(function(i) {rootB$node$AddDatum(i)
            self$assignments[[i]] <- rootB$node}, dataIdsA)
          
          # swap mains
          if (swapMain) {
            Descend <- function(root, nodeNum) {
              if (nodeNum == nodeAnum) root$main <- mainB
              if (nodeNum == nodeBnum) root$main <- mainA
              actChildren <- root$children[root$keep]
              for (i in seq_along(actChildren)) {
                nodeNum <- nodeNum + 1
                res <- Descend(actChildren[[i]], nodeNum)
                root$children[root$keep][[i]] <- res$root
                nodeNum <- res$nodeNum
              }
              return(list(root=root, nodeNum=nodeNum))
            }
            
            if (rootA$main != 1 && rootB$main != 1) {
              mainA <- rootA$main
              mainB <- rootB$main
              self$root <- Descend(root=self$root, nodeNum=1)$root
            }
          }
        }
        
        SwapNodes(nodeAnum, nodeBnum)
        self$ResampleSticks()
        post_new <- self$GetUnnormalizedPost()$ll
        if (post_new < post) {
          SwapNodes(nodeAnum, nodeBnum)
          self$ResampleSticks()
        } else {if (verbose) cat("Successful swap!\n")}
      }
    }
  )
)





