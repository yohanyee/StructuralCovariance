# Base object types
as.RSC.object <- function(x, ...) {
  class(x) <- "RSC.object"
  attributes(x) <- c(attributes(x), 
                     list(created=date(), 
                          ...))
  return(x)
}

print.RSC.object <- function(x, ...) {
  cat("RSC object of no set type\n")
  cat("Attributes:\n")
  cat(paste(names(attributes(x)), collapse=", "))
  cat("\n")
}

# Matrix types
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

as.RSC.cmatrix <- function(x, ...) {
  x <- as.RSC.matrix(x, RSC.matrix.type="correlation", ...)
  class(x) <- c("RSC.cmatrix", class(x))
  return(x)
}

as.RSC.pmatrix <- function(x, ...) {
  x <- as.RSC.matrix(x, 
                     RSC.matrix.type="p.value", 
                     p.adjust.method="none", ...)
  class(x) <- c("RSC.pmatrix", class(x))
  return(x)
}

print.RSC.matrix <- function(x, ...) {
  switch(as.character(attr(x, "RSC.matrix.type")),
         correlation={ cat("RSC correlation matrix\n") },
         p.value={ cat("RSC p.value matrix\n") }, 
         `NA`={ cat("RSC matrix of no set type\n") },
         { cat("RSC matrix of undefined type\n") }
         )
  cat(paste("Size:", paste(attr(x, "RSC.matrix.size"), collapse=" x "), "\n"))
  cat(paste("Ordered:", attr(x, "RSC.matrix.ordered"), "\n"))
}

# List types 
as.RSC.list <- function(x, ...) {
  x <- as.list(x)
  x <- as.RSC.object(x, 
                     RSC.object.type="list", 
                     RSC.list.type=NA,
                     RSC.list.size=length(x), ...)
  class(x) <- c("RSC.list", class(x))
  return(x)
}

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

as.RSC.speciallist <- function(x, ...) {
  x <- as.RSC.list(x, RSC.list.type="special", ...)
  class(x) <- c("RSC.speciallist", class(x))
  return(x)
}

print.RSC.list <- function(x, ...) {
  switch(as.character(attr(x, "RSC.list.type")),
         matrix={ cat("RSC matrix list\n") }, 
         `NA`={ cat("RSC list of no set type\n") },
         { cat("RSC list of undefined type\n") }
  )
  cat(paste("Size:", attr(x, "RSC.list.size"), "\n"))
  cat("Contents:\n")
  if (is.null(attr(x, "names"))) {
    cat(paste(attr(x, "RSC.list.size"), "unnamed items\n"))
  } else {
    if (attr(x, "RSC.list.size") >= 10) {
      cat(paste(attr(x, "names")[1:5], collapse=", "))
      cat(" ... ")
      cat(paste("<", as.numeric(attr(x, "RSC.list.size"))-6, "more >"))
      cat(" ... ")
      cat(paste(rev(attr(x, "names"))[1], "\n"))
    } else {
      cat(paste(attr(x, "names"), collapse=", "))
      cat("\n")
    }
  }
}

# TODO
#
# Special types
# 
# as.RSC.special
# 
# as.RSC.cormatrix
# 
# as.RSC.multicormatrix
# 
# print.RSC.special



