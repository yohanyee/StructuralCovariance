# Base object types

#' @export
as.RSC.object <- function(x, ...) {
  class(x) <- "RSC.object"
  attributes(x) <- c(attributes(x), 
                     list(created=date(), 
                          ...))
  return(x)
}

# Matrix types

#' @export
as.RSC.matrix <- function(x, ...) {
  x <- as.matrix(x)
  x <- as.RSC.object(x, 
                     RSC.object.type="matrix", 
                     RSC.matrix.type=NA,
                     RSC.matrix.size=attr(x,"dim"),
                     RSC.matrix.ordered=FALSE, ...)
  class(x) <- c("RSC.matrix", class(x))
  return(x)
}

#' @export
as.RSC.cmatrix <- function(x, ...) {
  x <- as.RSC.matrix(x, RSC.matrix.type="correlation", ...)
  class(x) <- c("RSC.cmatrix", class(x))
  return(x)
}

#' @export
as.RSC.pmatrix <- function(x, ...) {
  x <- as.RSC.matrix(x, 
                     RSC.matrix.type="p.value", 
                     p.adjust.method="none", ...)
  class(x) <- c("RSC.pmatrix", class(x))
  return(x)
}

# List types 

#' @export
as.RSC.list <- function(x, ...) {
  x <- as.list(x)
  x <- as.RSC.object(x, 
                     RSC.object.type="list", 
                     RSC.list.type=NA,
                     RSC.list.size=length(x), ...)
  class(x) <- c("RSC.list", class(x))
  return(x)
}

#' @export
as.RSC.matrixlist <- function(x, ...) {
  x <- as.RSC.list(x, RSC.list.type="matrix", ...)
  if (!all(sapply(x, attr, "RSC.object.type")=="matrix")) {
    tryCatch({x <- lapply(x, as.RSC.matrix)
    warning("One or more matrices did not have RSC.object.type of \"matrix\"; coerced all matrices to RSC.matrix")
    }, error=function(e) {
      stop("All matrices in list must have RSC.object.type of \"matrix\"")         
    })
  }
  class(x) <- c("RSC.matrixlist", class(x))
  return(x)
}

#' @export
as.RSC.speciallist <- function(x, ...) {
  x <- as.RSC.list(x, RSC.list.type="special", ...)
  class(x) <- c("RSC.speciallist", class(x))
  return(x)
}


# Special types

#' @export
as.RSC.special <- function(x, ...) {
  x <- as.RSC.object(x, 
                     RSC.object.type="special", 
                     RSC.special.type=NA, ...)
  class(x) <- c("RSC.special", class(x))
  return(x)
}

#' @export
as.RSC.cormatrix <- function(x, ...) {
  if (!all(c("correlation", "p.value") %in% names(x))) {
    stop("List must contain correlation and p.value matrices")
  }
  if (class(x$correlation)[1]!="RSC.cmatrix" | attr(x$correlation, "RSC.matrix.type") != "correlation") {
    stop("Correlation matrix is not the right type")
  }
  if (class(x$p.value)[1]!="RSC.pmatrix" | attr(x$p.value, "RSC.matrix.type") != "p.value") {
    stop("P.value matrix is not the right type")
  }
  x <- as.RSC.matrixlist(list(correlation=x$correlation, p.value=x$p.value), 
                     RSC.object.type="special", 
                     RSC.special.type="cormatrix", ...)
  class(x) <- c("RSC.cormatrix", class(x))
  return(x)
}

# TODO
#
# as.RSC.multicormatrix
# 
# print.RSC.special



