# MRSL: A phenome-wide causal discovery algorithm based on summarized genetic data

The marginal causal graph can be ontained by STEP1(), default method is IVW, change the method as MREgger by setting use_eggers_step1=T
STEP1(data_sum_beta,data_sum_se,beta,use_eggers_step1=F,p_cut=0.05)

The conditional causal graph can be obtained by STEP2(), default method is IVW, change the method as MREgger by setting use_eggers_step1=T;
STEP2(amatt,data_sum_beta,data_sum_se,beta_new,adj_methods=1,use_eggers_step2=F,vary_mvmr_adj=0,coefi=coeft_new,p_cut=0.05)
adj_methods=1,adjusting for all the nodes except colloders;adj_methods=2,adjusting for all the nodes on the open paths;adj_methods=3,adjusting for the minimal sufficient adjustment set and all the mediators.
