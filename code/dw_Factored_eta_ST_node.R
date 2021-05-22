#' Node is a R6 object of each cluster in TSSB
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords Node class
#' @field dataIds: data IDs
#' @field tssb: A TSSB object
#' @field children: descentant node objects
#' @field parent: parent node object
#' @method new
#' @method GetChildren
#' @method GetParent
#' @method SetParent
#' @method AddChild
#' @method RemoveChild
#' @method Kill
#' @method Spawn
#' @method HasData
#' @method AddDatum
#' @method RemoveDatum
#' @method GetNumOfLocalData
#' @method GetNumOfSubTreeData
#' @method GetData
#' @method GetAncestors

Node_DW_Factored_eta_ST <- R6::R6Class(
  classname = "Node_DW_Factored_eta_ST",

  public = list(
    # Fields ------------------------------------------------------------------
    dataIds = c(),
    tssb = NULL,
    numOfChildren = NULL,  # New added!!
    dataDims = NULL,
    depth = NULL,

    # Methods -----------------------------------------------------------------
    initialize = function(parent = NULL, 
                          tssb = NULL, 
                          numOfChildren = NULL,
                          dataDims = NULL,
                          depth = NULL) {
      if (!is.null(parent)) {
        private$parent <- parent
        parent$AddChild(self)
      } else {
        private$parent <- parent
      }
      if (!is.null(tssb) ) self$tssb <- tssb
      
      # New added!
      self$numOfChildren <- numOfChildren 
      self$dataDims <- dataDims
      self$depth <- depth
      },


    # 因为 children 是私有成员，所以外部不能 node$children 的方式调用
    # 但 GetChildren 是公有方法，所以可以通过 node$GetChildren 的方式取出 children！！
    # 疑问：为何非要把 children 和 parents 设置为私有成员呢？？？
    GetChildren = function() {
      private$children
    },

    GetParent = function() {
      private$parent
    },

    SetParent = function(parent = NULL) {
      private$parent <- parent
    },

    AddChild = function(child) {
      if (!missing(child)) {
        child$SetParent(self)
        private$children <- c(private$children, child)
      }
      invisible(self)
    },

    RemoveChild = function(child) {
      if (!missing(child)) {
        private$children <- unlist(Filter(
          Negate(  
            # Negate：取反
            # 此处是取出和 child 不同的所有 children，重新赋值给 private，达到 Remove child 的目的
            function(x) identical(child, x)
            ),
          private$children))
      }
      invisible(self)
    },

    Kill = function() {
      if (!is.null(private$parent))  {
      # 不会 kill 掉 root node！
        private$parent$RemoveChild(self) # 它的父亲把作为孩子的它杀掉了
      }
      private$children = NULL  # 被父亲杀掉之后，杀掉自己的孩子
      private$parent = NULL    # 断绝自己的父亲
    },

    Spawn = function() {
      # 当产生新的 child node 时：
      # parent：为实例化的 self
      # tssb：为实例化 self 所在的 tssb
      # numOfChildren：和 self 相同，都是 maxWidth
      return(get(class(self)[1])$new(parent = self, tssb = self$tssb, 
                                     numOfChildren = self$numOfChildren,
                                     dataDims = self$dataDims,
                                     depth = self$depth+1))
      # get 函数：Return the Value of a Named Object
      # 此处返回：<Normal_DW> object generator（生成器）
    },

    # 判断是否有数据位于此 node 或经过此 node
    HasData = function() {
      if ( length(self$dataIds)>0 ) {
        return(TRUE)
      } else {
        return(
          sum(
            unlist(sapply(private$children, function(x) x$HasData()))
            )
          > 0)
      }
    },

    AddDatum = function(id) {
      if (!missing(id) && !any(self$dataIds==id)) {
        self$dataIds <- c(self$dataIds, id)
      }
    },

    RemoveDatum = function(id) {
      if (!missing(id)) {
        if (!id %in% self$dataIds) {
          warning("id is not found in dataIds, nothing is removed")
        } else {
          self$dataIds <- self$dataIds[self$dataIds != id]
        }
      }
    },

    GetNumOfLocalData = function() {
      length(self$dataIds)
    },

    GetNumOfSubTreeData = function() {
        Reduce(
          sum,
          Map(function(x) x$GetNumOfSubTreeData(),
              private$children),
          length(self$dataIds)
        )
    },

    GetData = function() {
      if (dim(self$tssb$data)[2] > 1) {
        return(self$tssb$data[self$dataIds,])
      } else {
        return(as.matrix(self$tssb$data[self$dataIds,]))
      }

    },

    GetAncestors = function() {
      ancestors = c()
      if (is.null(private$parent)) {
        return(list(self))
      } else {
        ancestors <- c(private$parent$GetAncestors(), self)
        return(ancestors)
      }
    }

    ), # end of public

  private = list(
    # Fields
    # 因为是私有成员，因此在 tssb 类中无法查看 node 的这两个成员！
    children = NULL,
    parent = NULL
    )
  )




