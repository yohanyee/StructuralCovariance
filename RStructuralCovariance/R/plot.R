#' @export
getColorFunction <- function(x, n=100) {
  UseMethod("getColorFunction", x)
}

#' @export
getColorFunction.default <- function(x, n) {
  cmap <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
  return(cmap(n))
}

#' @export
getColorFunction.RSC.cmatrix <- function(x, n) {
  cmap <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
  return(cmap(n))
}

#' @export
getColorFunction.RSC.pmatrix <- function(x, n) {
  cmap <- colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  return(cmap(n))
}

#' @export
plot.RSC.matrix <- function(mat, package="corrplot", ...) {
  switch(package,
         corrplot={
           plt <- corrplot::corrplot(mat, method="color", col=getColorFunction(mat), ...)
         },
         gplots={
           plt <- gplots::heatmap.2(mat, 
                                    Rowv=FALSE, colV=FALSE,
                                    trace="none",
                                    ...
           )
         },
         ggplot={cat("Not implemented yet.\n")}
         )
  print(plt)
}

#plot.RSC.matrixlist

#plot.RSC.cormatrix