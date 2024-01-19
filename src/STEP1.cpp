#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include <R.h>
#include <iostream>

using namespace arma;
using namespace std;
using namespace Rcpp;


// combination
bool next_comb(int* comb, const int n, const int k) {
  int i = k - 1;
  const int e = n - k;
  do
    comb[i]++;
  while (comb[i] > e + i && i--);
  if (comb[0] > e)
    return 0;
  while (++i < k)
    comb[i] = comb[i - 1] + 1;
  return 1;
}

int Combination(int* cho, int nl, int n, int k) {

  if (n < k || k <= 0)
    return 0;
  int* comb = new int[k];
  for (int i = 0; i < k; i++)
    comb[i] = i;

  do
    for (int i = 0; i < k; ++i) {
      cho[nl] = comb[i];
      //cout << cho[nl] << endl;
      nl++;
    }
    while (next_comb(comb, n, k));

  delete[] comb;
  return nl;
}


// choose IVs
void Cho_IV(const int nodex, const arma::mat& bet,arma::uvec& l1){

  vec ttsnp(bet.n_rows);
  for(int rsnp=0; rsnp<bet.n_rows; rsnp++){
    if(sum(bet.row(rsnp)!=0)==1){
      ttsnp[rsnp]=1;
    }else{
      ttsnp[rsnp]=0;
    }
  }

  mat bet0=bet.rows(find(ttsnp == 1));
  l1=find(bet0.col(nodex)!= 0);

 // cout << "l1==" << l1 << endl;

}


// linear regression
Rcpp::List fastLm(const arma::mat& X1, const arma::colvec& y,bool intercept) {

  int n = X1.n_rows;

  mat X;
  if(intercept)
    X=join_rows(ones(n,1),X1);
  else
    X=X1;

  int k = X.n_cols;

  //arma::colvec coef = arma::spsolve(X, y);    // fit model y ~ X
  arma::mat coef = arma::pinv(X) * y;
  arma::colvec res  = y - X*coef;           // residuals

  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);

  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

  arma::colvec tval = coef/std_err;

  NumericVector xx=as<NumericVector>(wrap(-abs(tval)));
  double df=n-k;

  NumericVector pval = Rcpp::pt(xx,df)*2;

  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr")       = std_err,
                            Rcpp::Named("tvalue")       = tval,
                            Rcpp::Named("pvalue")       = pval,
                            Rcpp::Named("df")  = n - k);
}


// judge edge
int Judge(const arma::mat& betaX,const arma::vec& betaY,const arma::vec& seY,
          const int use_egger){

  int edge=0;
  mat sigma = diagmat(1/seY);
  mat betaX1 = sigma*betaX;
  mat betaY1 = sigma*betaY;

  if(use_egger==1){

    List lm1=fastLm(betaX1,betaY1,true);
    arma::vec pval1 = as<arma::vec>(lm1["pvalue"]);

    if(pval1[0]>0.05){
      List lm2=fastLm(betaX1,betaY1,false);
      arma::vec pval2 = as<arma::vec>(lm2["pvalue"]);
      if(pval2[0]<0.05)
        edge=1;
    }else if(pval1[1]<0.05){
      edge=1;
    }

  }else if(use_egger==0){

    List lm2=fastLm(betaX1,betaY1,false);
    arma::vec pval2 = as<arma::vec>(lm2["pvalue"]);
    if(pval2[0]<0.05)
      edge=1;

  }

  return edge;

}


// [[Rcpp::export]]

SEXP STEP1(SEXP data_sum_beta,SEXP data_sum_se,SEXP beta1,
            SEXP use_eggers_step1){

  const arma::mat sum_beta = as<arma::mat>(data_sum_beta);
  const arma::mat sum_se = as<arma::mat>(data_sum_se);
  const arma::mat beta = as<arma::mat>(beta1);
  const int use_egger_step1=as<int>(use_eggers_step1);

  int D1 = sum_beta.n_cols, D2 = sum_se.n_cols;
  if (D1 != D2){
    perror("The dimensions of beta and se are not matched");
  }

  mat amatt;
  amatt = zeros<mat>(D1,D1);

  int* cho = new int[100000000];
  int nl=0;
  nl=Combination(cho,nl,D1,2);

  // step 1
  cout << "step1 all for the first direction " << nl << "..." << endl;
  for(int i=0; i<nl/2; i++){
    // int i=0;
    uvec l1;
    Cho_IV(cho[2*i],beta,l1);

    //cout << cho[2*i] << " " << cho[2*i+1] << endl;
    mat betaX1 = sum_beta.col(cho[2*i]);
    betaX1 = betaX1.elem(l1);
    vec betaY1 = sum_beta.col(cho[2*i + 1]);
    betaY1 = betaY1.elem(l1);
    vec seY1 = sum_se.col(cho[2*i + 1]);
    seY1 = seY1.elem(l1);
    amatt(cho[2*i],cho[2*i + 1]) = Judge(betaX1,betaY1,seY1,use_egger_step1);

  }

  cout << "step1 all for the second direction " << nl << "..." << endl;
  for(int i=0; i<nl/2; i++){
    // int i=0;
    uvec l2;
    Cho_IV(cho[2*i+1],beta,l2);

    //cout << cho[2*i+1] << " " << cho[2*i] << endl;
    mat betaX2 = sum_beta.col(cho[2*i+1]);
    betaX2 = betaX2.elem(l2);
    vec betaY2 = sum_beta.col(cho[2*i]);
    betaY2 = betaY2.elem(l2);
    vec seY2 = sum_se.col(cho[2*i]);
    seY2 = seY2.elem(l2);
    amatt(cho[2*i+1],cho[2*i]) = Judge(betaX2,betaY2,seY2,use_egger_step1);

  }

  cout << "step1 over" << endl;
  //cout << "amatt " << amatt << endl;

  //cout << "memory " << GetSysMemInfo() << endl;


  delete[] cho;

  return List::create(Rcpp::Named("amatt") = amatt);
}
