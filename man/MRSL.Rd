% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRSL.R
\name{MRSL}
\alias{MRSL}
\title{MRSL}
\usage{
MRSL(
  M_amatt,
  data_sum_beta,
  data_sum_se,
  beta,
  cutoff,
  adj_methods,
  use_eggers_step2,
  vary_mvmr_adj
)
}
\arguments{
\item{M_amatt}{The D*D adjacent matrix of marginal causal graph.}

\item{data_sum_beta}{J*D matrix for summary coefficients of all J SNPs
on the D variables in the marginal causal graph, e.g. coefficients
from linear regression of each variable on each SNP.}

\item{data_sum_se}{J*D matrix of standard error corresponding to
data_sum_beta.}

\item{beta}{J*D matrix of indicator with element 0/1,
e.g. the element in the i-th row and the j-th column is 1 denote
the i-th SNP is the IV for the j-th node.}

\item{cutoff}{The cutoff value of p value to determine which
edges' effect is significant.}

\item{adj_methods}{Three methods for choosing adjusting variables:
1 remove collier; 2 variables in open path; 3 min sufficient sets}

\item{use_eggers_step2}{MVMR methods used in the second step:
1 MVMR-Egger; 0 MVMR-IVW}

\item{vary_mvmr_adj}{vary_mvmr_adj=0 denotes directly adjusting for
the adjusting variables in MVMR;vary_mvmr_adj=1 denotes varying
the number of adjusting variables in MVMR from the selected
variable by adj_methods.}
}
\value{
A list object that contains a summary of conditional causal graph,
    amatt is the adjacent matrix of conditional causal graph;
    iteration is the number of iteration in the third step of MRSL;
    nnames is the topological sorting of nodes.
}
\description{
Provides an adjacent matrix of conditional causal graph
    Elements of adjacent matrix are 0/1, 1 denotes the edge is exist.
}
\examples{
load('toyexample.Rdata')
res <- MRSL(M_amatt,data_sum_beta,data_sum_se,beta,cutoff=0.05,
             adj_methods=1,use_eggers_step2=0,vary_mvmr_adj=0)
res$amatt

}
