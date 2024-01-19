#' @title PLOT_NET
#'
#' @description Plot a netwok using an adjacent matrix of a graph
#'
#' @param amatt An adjacent matrix of a graph
#'
#' @return A value of Kendall’s tau
#'
#' @examples
#' gg <- PLOT_NET(res$amatt)
#' V(gg)$size = 18
#' plot(gg,layout=layout_with_kk(gg))
#' # the output of PLOT_NET can be operated as the igraph
#'
#' @export
#'
PLOT_NET <- function(amatt){

  nnames <- colnames(amatt)
  ggrap <- nnames[c(t(which(amatt!=0,arr.ind = T)))]
  gg <- graph(ggrap)

  return(gg)
}
