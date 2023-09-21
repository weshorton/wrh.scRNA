comparePCs <- function(obj1, obj2, dr1 = "pca", dr2 = NULL, dim1, dim2 = NULL, nFeat_v = 100, which_v = "all") {
  #' Compare top n genes of PCs from different seurat objects
  #' @description 
  #' Compare the top genes of a given PC from two different seurat objects
  #' @param obj1 First seurat object
  #' @param obj2 Second seurat object
  #' @param dr1 Name of dimensional reduction to grab for obj1. Default is 'pca'
  #' @param dr2 Name of dimensional reduction to grab fro obj2. If NULL (default) will be same as dr1.
  #' @param dim1 Index of principal component to use for obj1
  #' @param dim2 Index of principal component to use for obj2. If NULL (default) will be same as dim1.
  #' @param nFeat_v Number of features to compare.
  #' @param which_v Which genes to compare. "all" takes both positive and negative genes; 
  #' "positive" is positive only; "negative" is negative only
  #' @return Vector of gene names shared between the two PCs
  #' @export
  
  ### Handle null arguments
  if (is.null(dim2)) dim2 <- dim1
  if (is.null(dr2)) dr2 <- dr1
  
  ### Handle which_v
  if (which_v == "all") {
    balance_v <- F
  } else {
    balance_v <- T
    nFeat_v <- 2*nFeat_v
  } # fi
  
  ### Get Features
  feat1 <- Seurat::TopFeatures(object = obj1[[dr1]],
                               dim = dim1,
                               nfeatures = nFeat_v,
                               balanced = balance_v)
  
  feat2 <- Seurat::TopFeatures(object = obj2[[dr2]],
                               dim = dim2,
                               nfeatures = nFeat_v,
                               balanced = balance_v)
  
  ### Get either positive, negative, or both.
  if (which_v == "positive") {
    feat1 <- feat1$positive
    feat2 <- feat2$positive
  } else if (which_v == "negative") {
    feat1 <- feat1$negative
    feat2 <- feat2$negative
  }  # fi
  
  ### Compare
  shared_v <- intersect(feat1, feat2)
  
  ### Output
  return(shared_v)
  
} # comparePCs
