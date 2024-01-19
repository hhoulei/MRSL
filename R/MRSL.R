#' @title MRSL
#'
#' @description Provides an adjacent matrix of conditional causal graph
#'     Elements of adjacent matrix are 0/1, 1 denotes the edge is exist.
#'
#' @param M_amatt The D*D adjacent matrix of marginal causal graph.
#' @param data_sum_beta J*D matrix for summary coefficients of all J SNPs
#'      on the D variables in the marginal causal graph, e.g. coefficients
#'      from linear regression of each variable on each SNP.
#' @param data_sum_se J*D matrix of standard error corresponding to
#'      data_sum_beta.
#' @param beta J*D matrix of indicator with element 0/1,
#'      e.g. the element in the i-th row and the j-th column is 1 denote
#'      the i-th SNP is the IV for the j-th node.
#' @param cutoff The cutoff value of p value to determine which
#'      edges' effect is significant.
#' @param adj_methods Three methods for choosing adjusting variables:
#'      1 remove collier; 2 variables in open path; 3 min sufficient sets
#' @param use_eggers_step2 MVMR methods used in the second step:
#'      1 MVMR-Egger; 0 MVMR-IVW
#' @param vary_mvmr_adj vary_mvmr_adj=0 denotes directly adjusting for
#'      the adjusting variables in MVMR;vary_mvmr_adj=1 denotes varying
#'      the number of adjusting variables in MVMR from the selected
#'      variable by adj_methods.
#' @useDynLib MRSL, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @return A list object that contains a summary of conditional causal graph,
#'     amatt is the adjacent matrix of conditional causal graph;
#'     iteration is the number of iteration in the third step of MRSL;
#'     nnames is the topological sorting of nodes.
#'
#'
#' @examples
#' load('toyexample.Rdata')
#' res <- MRSL(M_amatt,data_sum_beta,data_sum_se,beta,cutoff=0.05,
#'              adj_methods=1,use_eggers_step2=0,vary_mvmr_adj=0)
#' res$amatt
#'
#' @export
#'
#'
#'
MRSL <- function(M_amatt,data_sum_beta,data_sum_se,beta,cutoff,
                 adj_methods,use_eggers_step2,vary_mvmr_adj){

  amatt0 <- list()
  amatt0[[1]] <- M_amatt

  n_nodes <- ncol(M_amatt)
  nnames <- colnames(M_amatt)

  # topological sorting
  pp <- amatt0[[1]]
  v1 <- n_nodes-c(DFS(amattt=pp,n_nodes=n_nodes))

  amatt0[[1]] <- amatt0[[1]][v1,v1]
  nnames <- nnames[v1]

  # iteration of MVMR
  just <- TRUE
  i=1
  while(just){
    cat("iteration--",i)
    i=i+1
    amatt055 <- STEP2(amatt0[[i-1]],data_sum_beta,data_sum_se,beta,
                      adj_methods,use_eggers_step2,vary_mvmr_adj)
    amatt0[[i]] <- amatt055$amatt
    just <- !all(amatt0[[i]]==amatt0[[i-1]])
  }
  amatt=amatt0[[i]]
  colnames(amatt) <- rownames(amatt) <- nnames

  return(list(amatt=amatt,
              iteration=i,
              nnames=nnames))
}
