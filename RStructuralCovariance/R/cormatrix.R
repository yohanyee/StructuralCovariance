#' Construct a correlation matrix by correlating volumes
#' 
#' Creates a correlation matrix by computing pairwise correlations between volumes (column vectors of a \code{data.frame})
#' 
#' @param v volume \code{data.frame} or \code{RMINC::anatGetAll} output consisting of individuals as rows and volumes as columns
#' @param normalize logical value denoting whether volumes should be normalized (rows summed to 1)
#' @param cor.method character string giving the type of correlation coefficient being computed
#' @param p.adjust.method character string giving method to use when correcting p values for multiple comparisons
#' @param resample logical value denoting whether to resample individuals from the volume \code{data.frame}
#' @param resample_size integer value between 3 and the number of rows of the volume \code{data.frame} 
#' @param resample_replace logical value denoting whether to resample with replacement
#' @param quiet logical value denoting whether progress should be printed
#' 
#' @return an \code{RSC.cormatrix} object
#' @export
construct_cormatrix <- function(v, normalize=TRUE, 
                                cor.method="pearson", cor.alternative_hypothesis="two.sided", p.adjust.method="fdr", 
                                resample=FALSE, resample_size=NULL, resample_replace=TRUE,
                                quiet=FALSE, ...) {
  v <- (if (normalize) { normalize_volumes(v) } else v)
  v <- (if (resample) { resample_volumes(v, size=resample_size, replace=resample_replace) } else v)
  cormatrix <- correlate_volumes(v, method=cor.method, alternate=cor.alternate_hypothesis, quiet=quiet)
  cormatrix <- p.adjust(cormatrix, method=p.adjust.method)
  return(cormatrix)
}

