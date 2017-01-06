#' @export
print.RSC.object <- function(x, ...) {
  cat("RSC object of no set type\n")
  cat("Attributes:\n")
  cat(paste(names(attributes(x)), collapse=", "))
  cat("\n")
}

#' @export
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

#' @export
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