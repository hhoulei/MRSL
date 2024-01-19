# MRSL: A Causal Network Pruning Algorithm based on GWAS Summary Data
 
The workflow of this algorithm is as below (Confounders of d phenotypes are omitted.).  

![Figure 1](https://github.com/hhoulei/MRSL_simul/blob/main/Figure%201_00.png)  

The input includes GWAS summary data for each phenotype and marginal causal graph. For the Step 1, the topological sorting of marginal causal graph should be found using Depth First Search (DFS). For the Step 2, MVMR is performed to remove extra edges in the marginal causal graph by adjusting for the genetic associations with phenotypes of three strategies of sufficient separating sets. Then a conditional causal graph is obtained. The in the Step 3, iteration for step 2 is performed until the graph converges. Finally, MRSL outputs an estimated conditional causal graph. (A-J) Motivating example with five nodes. (A) The true causal graph. (B) Marginal causal graph. (C-J) Perform MVMR for each edge in graph (F) based on its topological sorting. Blue nodes denote the exposure and outcome we are interested in. MRSL outputs the graph (J). 

***Parameters Input***   
`M_amatt` The D×D adjacent matrix of marginal causal graph.  
`data_sum_beta` J×D matrix for summary coefficients of all J SNPs on the D variables in the marginal causal graph, e.g. coefficients from linear regression of each variable on each SNP.  
`data_sum_se` J×D matrix of standard errors corresponding to data_sum_beta.  
`beta` J×D matrix of indicators with element 0/1, e.g. the element in the i-th row and the j-th column is 1 denote the i-th SNP is the IV for the j-th node.  
`cutoff` The cutoff value of p-value determining which edges' effect is significant.  
`adj_methods` Three methods for choosing adjusting variables: 1 remove collier; 2 variables in open path; 3 min sufficient sets  
`use_eggers_step2` MVMR methods used in the second step: 1 MVMR-Egger; 0 MVMR-IVW  
`vary_mvmr_adj` vary_mvmr_adj=0 denotes directly adjusting for the adjusting variables in MVMR; vary_mvmr_adj=1 denotes varying the number of adjusting variables in MVMR from the selected variable by adj_methods.  

***Output***   
A list object that contains a summary of conditional causal graph,   
`amatt` is the adjacent matrix of conditional causal graph;  
`iteration` is the number of iteration in the third step of MRSL;  
`nnames` is the topological sorting of nodes.  

***Installation***  
`devtools::install_github("hhoulei/MRSL")`  
This installation depends on the compilation environment in your systerm (LAPACK and BLAS libraries). Possible compilation problems and workarounds can be found in the Issues section (https://github.com/hhoulei/MRSL/issues). We welcome you to continue to discuss your errors and solutions.  

***Toy Example***  
The following toyexample.Rdata can be downloaded from https://github.com/hhoulei/MRSL_simul/blob/main/toyexample.Rdata   
`library(MRSL)`  
`library(igraph)`  
`load('toyexample.Rdata')`  
`res <- MRSL(M_amatt,data_sum_beta,data_sum_se,`  
`beta,cutoff=0.05,adj_methods=1,use_eggers_step2=0,vary_mvmr_adj=0)`  
`res$amatt`  
plot the conditional causal graph use PLOT_NET:  
`gg <- PLOT_NET(res$amatt)`  
`plot(gg,layout=layout_with_kk(gg))`  


***Citation***:  
MRSL: A phenome-wide causal discovery algorithm based on GWAS summary data  
Lei Hou, Zhi Geng, Xu Shi, Chuan Wang, Hongkai Li, Fuzhong Xue  
medRxiv 2022.06.29.22277051; doi: https://doi.org/10.1101/2022.06.29.22277051  

Please contact houlei@pku.edu.cn for any questions. We will continue to update this R package and reduce the problems that may be encountered during its installation.
