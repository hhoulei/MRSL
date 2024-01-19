#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include <R.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;


arma::vec Open_Path_Adj(const mat& MM, const int start, const int end,
                        const int adj){

  Function f("Open_path_adj");
  NumericVector loc_adj =  f(Named("MM")=MM,
                    Named("start")=start,
                    Named("end")=end,
                    Named("adj")=adj);

  arma::vec loc_adj1 = as<arma::vec>(loc_adj);

  //cout << "loc_adj1==" << loc_adj1 << endl;

  return loc_adj1;
}


// combination
bool next_comb2(int* comb, const int n, const int k) {
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

int Combination2(int* cho, int nl, int n, int k) {

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
    while (next_comb2(comb, n, k));

  delete[] comb;
  return nl;
}


// choose IVs
void Cho_IV1(const int nodex1, const arma::vec& nodeadjust, const arma::mat& bet1,arma::uvec& l19){

  uvec posadj(nodeadjust.n_elem+1);
  for(int hu=0; hu<(nodeadjust.n_elem); hu++){
    posadj[hu]=nodeadjust[hu];
  }
  posadj[nodeadjust.n_elem]=nodex1;

  mat mstadj = bet1.cols(posadj);
  vec ll1(mstadj.n_rows);
  for(int ko=0; ko<(mstadj.n_rows); ko++){
    if(any(mstadj.row(ko)!= 0)) ll1[ko]=1;
  }
  l19=find(ll1==1);

}


// linear regression
Rcpp::List fastLm2(const arma::mat& X1, const arma::colvec& y,bool intercept) {

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
int Judge2(const arma::mat& betaX,const arma::vec& betaY,const arma::vec& seY,
          const int use_egger){

  int edge=0;
  mat sigma = diagmat(1/seY);
  mat betaX1 = sigma*betaX;
  mat betaY1 = sigma*betaY;

  if(use_egger==1){

    List lm1=fastLm2(betaX1,betaY1,true);
    arma::vec pval1 = as<arma::vec>(lm1["pvalue"]);

    if(pval1[0]>0.05){
      List lm2=fastLm2(betaX1,betaY1,false);
      arma::vec pval2 = as<arma::vec>(lm2["pvalue"]);
      if(pval2[0]<0.05)
        edge=1;
    }else if(pval1[1]<0.05){
      edge=1;
    }

  }else if(use_egger==0){

    List lm2=fastLm2(betaX1,betaY1,false);
    arma::vec pval2 = as<arma::vec>(lm2["pvalue"]);
    if(pval2[0]<0.05)
      edge=1;

  }

  return edge;

}


// [[Rcpp::export]]

SEXP STEP2(SEXP amattt, SEXP data_sum_beta,SEXP data_sum_se,SEXP beta1,
           SEXP adj_methods,SEXP use_eggers_step2,SEXP vary_mvmr_adj){

  arma::mat amatt = as<arma::mat>(amattt);
  const arma::mat sum_beta = as<arma::mat>(data_sum_beta);
  const arma::mat sum_se = as<arma::mat>(data_sum_se);
  const arma::mat beta = as<arma::mat>(beta1);

  const int adj_method=as<int>(adj_methods);
  const int use_egger_step2=as<int>(use_eggers_step2);
  const int vary_mvmr_adjnodes=as<int>(vary_mvmr_adj);


  int D1 = sum_beta.n_cols, D2 = sum_se.n_cols;
  if (D1 != D2){
    perror("The dimensions of beta and se are not matched");
  }


  int* cho = new int[100000000];
  int nl=0;
  nl=Combination2(cho,nl,D1,2);

  // step 2

  cout << "step2 all start the first direction--" << nl << "..." << endl;

  for(int i=0; i<(nl/2); i++){

    //cout << "step2--first " << i << endl;

    arma::vec nodesadj_all;

    if(amatt(cho[2*i],cho[2*i + 1])==1){

      if(adj_method==1){

        vec allnodes=regspace(0,D1-1);
        // cout << allnodes.n_elem << endl;

        // remove itself
        uvec remit={0,0};;
        remit[0]=cho[2*i];
        remit[1]=cho[2*i+1];
        allnodes.shed_rows(remit);

        //cout << "remove colliders" << endl;
        vec rem(D1-2);
        for(int k=0; k<(D1-2); k++){
          rem[k]=amatt(cho[2*i],allnodes[k])+amatt(cho[2*i+1],allnodes[k]);
        }
        nodesadj_all=allnodes.elem(find(rem<2));

      }else{

        //cout << "problem" << endl;
        nodesadj_all = Open_Path_Adj(amatt,cho[2*i]+1,cho[2*i + 1]+1,
                                           adj_method);
        //cout << "cho[2*i]+1==" << cho[2*i]+1 << endl;
        //cout << "cho[2*i + 1]+1==" << cho[2*i + 1]+1 << endl;
      }

      // cout << "cpp--num of nodesadj==" << nodesadj.n_elem << endl;

      if(nodesadj_all.n_elem==0){
        //cout << "No variable for MVMR" << endl;
      }else{

        // MVMR

        int edge_up;

        //cout << "MVMR... " << endl;

        if(vary_mvmr_adjnodes==1){

          int nodesnum=nodesadj_all.n_elem;

          for(int k=0; k<nodesnum; k++){

            // int k=0;

            int* cho2 = new int[100000000];
            int nl2=0;
            nl2=Combination2(cho2,nl2,nodesnum,k+1);

            //cout << "nodesnum k -- " << k << " nl2 " << nl2 << endl;

            // cout << cho[2*i] << " " << cho[2*i+1] << endl;


            for(long e=0; e<(nl2/(k+1)); e++){

              vec nodesadj=nodesadj_all.rows(cho2[(k+1)*e],cho2[(k+1)*e+k]);
              //cout << "nodesadj.n_elem " << nodesadj.n_elem << endl;
              uvec l3;

              Cho_IV1(cho[2*i],nodesadj,beta,l3);


              mat betaX3 = sum_beta.col(cho[2*i]);
              betaX3 = betaX3.elem(l3);
              vec betaY3 = sum_beta.col(cho[2*i+1]);
              betaY3 = betaY3.elem(l3);
              vec seY3 = sum_se.col(cho[2*i+1]);
              seY3 = seY3.elem(l3);

              vec choi;
              vec betaX4;
              mat betaX5=betaX3;
              for(int ut=0; ut<(k+1); ut++){

                choi=nodesadj_all.row(cho2[(k+1)*e+ut]);
                int choi1=choi[0];
                //cout << "choi " << choi << endl;

                betaX4 = sum_beta.col(choi1);
                betaX4 = betaX4.elem(l3);
                betaX5=join_rows(betaX5,betaX4);
              }

              edge_up=Judge2(betaX5,betaY3,seY3,use_egger_step2);
              if(edge_up==0){
                amatt(cho[2*i],cho[2*i+1])=0;
                k=nodesnum;
                break;
              }
              // cout <<  "edge_up " << edge_up << endl;
            }
            delete[] cho2;
          }

        }else{
          uvec l3;

          Cho_IV1(cho[2*i],nodesadj_all,beta,l3);



          mat betaX3 = sum_beta.col(cho[2*i]);
          betaX3 = betaX3.elem(l3);
          vec betaY3 = sum_beta.col(cho[2*i+1]);
          betaY3 = betaY3.elem(l3);
          vec seY3 = sum_se.col(cho[2*i+1]);
          seY3 = seY3.elem(l3);

          vec betaX4;
          mat betaX5=betaX3;
          for(int ut=0; ut<(nodesadj_all.n_elem); ut++){

            int choi1=nodesadj_all[ut];
            betaX4 = sum_beta.col(choi1);
            betaX4 = betaX4.elem(l3);
            betaX5=join_rows(betaX5,betaX4);
          }

          edge_up=Judge2(betaX5,betaY3,seY3,use_egger_step2);
          if(edge_up==0){
            amatt(cho[2*i],cho[2*i+1])=0;
          }
        }



      }

    }

  }

  cout << "step2 all start the second direction" << nl << "..." << endl;
  for(int i=0; i<(nl/2); i++){

    //cout << "step2--second " << i << endl;
    arma::vec nodesadj_all;

    if(amatt(cho[2*i+1],cho[2*i])==1){


      if(adj_method==1){

        vec allnodes=regspace(0,D1-1);
        // cout << allnodes.n_elem << endl;

        // remove itself
        uvec remit={0,0};;
        remit[0]=cho[2*i];
        remit[1]=cho[2*i+1];
        allnodes.shed_rows(remit);

        //cout << "remove colliders" << endl;
        vec rem(D1-2);
        for(int k=0; k<(D1-2); k++){
          rem[k]=amatt(cho[2*i],allnodes[k])+amatt(cho[2*i+1],allnodes[k]);
        }
        nodesadj_all=allnodes.elem(find(rem<2));

      }else{

        nodesadj_all = Open_Path_Adj(amatt,cho[2*i]+1,cho[2*i + 1]+1,
                                           adj_method);
      }

      //cout << "cpp--num of nodesadj==" << nodesadj.n_elem << endl;

      if(nodesadj_all.n_elem==0){
        //cout << "No variable for MVMR" << endl;
      }else{

        //cout << "MVMR..." << endl;
        int edge_up;

        if(vary_mvmr_adjnodes==1){

          int nodesnum=nodesadj_all.n_elem;

          for(int k=0; k<nodesnum; k++){

            // int k=0;

            int* cho2 = new int[100000000];
            int nl2=0;
            nl2=Combination2(cho2,nl2,nodesnum,k+1);

            //cout << "nodesnum k -- " << k << " nl2 " << nl2 << endl;

            // cout << cho[2*i] << " " << cho[2*i+1] << endl;


            for(long e=0; e<(nl2/(k+1)); e++){

              vec nodesadj=nodesadj_all.rows(cho2[(k+1)*e],cho2[(k+1)*e+k]);
              //cout << "nodesadj.n_elem " << nodesadj.n_elem << endl;
              uvec l3;

              Cho_IV1(cho[2*i],nodesadj,beta,l3);


              mat betaX3 = sum_beta.col(cho[2*i]);
              betaX3 = betaX3.elem(l3);
              vec betaY3 = sum_beta.col(cho[2*i+1]);
              betaY3 = betaY3.elem(l3);
              vec seY3 = sum_se.col(cho[2*i+1]);
              seY3 = seY3.elem(l3);

              vec choi;
              vec betaX4;
              mat betaX5=betaX3;
              for(int ut=0; ut<(k+1); ut++){

                choi=nodesadj_all.row(cho2[(k+1)*e+ut]);
                int choi1=choi[0];
                //cout << "choi " << choi << endl;

                betaX4 = sum_beta.col(choi1);
                betaX4 = betaX4.elem(l3);
                betaX5=join_rows(betaX5,betaX4);
              }

              edge_up=Judge2(betaX5,betaY3,seY3,use_egger_step2);
              if(edge_up==0){
                amatt(cho[2*i],cho[2*i+1])=0;
                k=nodesnum;
                break;
              }
              // cout <<  "edge_up " << edge_up << endl;
            }
            delete[] cho2;
          }

        }else{

          // MVMR

          uvec l3;

          Cho_IV1(cho[2*i],nodesadj_all,beta,l3);


          mat betaX3 = sum_beta.col(cho[2*i]);
          betaX3 = betaX3.elem(l3);
          vec betaY3 = sum_beta.col(cho[2*i+1]);
          betaY3 = betaY3.elem(l3);
          vec seY3 = sum_se.col(cho[2*i+1]);
          seY3 = seY3.elem(l3);

          vec betaX4;
          mat betaX5=betaX3;
          for(int ut=0; ut<(nodesadj_all.n_elem); ut++){

            int choi1=nodesadj_all[ut];
            betaX4 = sum_beta.col(choi1);
            betaX4 = betaX4.elem(l3);
            betaX5=join_rows(betaX5,betaX4);
          }

          edge_up=Judge2(betaX5,betaY3,seY3,use_egger_step2);
          if(edge_up==0){
            amatt(cho[2*i],cho[2*i+1])=0;
          }

        }



      }

    }

  }


  delete[] cho;

  return List::create(Rcpp::Named("amatt") = amatt);
}
