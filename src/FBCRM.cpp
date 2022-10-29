#include <time.h>
#include <Rmath.h>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <iostream>     // std::cout
#include <cmath>
#include <cfloat>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//////Get absolute value for a single observation/////////////
//[[Rcpp::export]]
double abs1(double a){
  if(a<0){
    return(-a);
  }else{
    return(a);
  }
}


//Gets minimum of two values
//[[Rcpp::export]]
double min1(double a, double b){
  double z=0;
  if(a>=b){
    z=b;
  }else{
    z=a;
  }
  
  return(z);
}


//Gets maximum of two values
//[[Rcpp::export]]
double max1(double a, double b){
  double z=0;
  if(a>=b){
    z=a;
  }else{
    z=b;
  }
  
  return(z);
}


//Get the exact position upto which all the other positional values will be clustered below in FBCRM
//[[Rcpp::export]]
int upper_bound(arma::vec zeta, //zeta vector
                int a //starting point in cpp, starts from 0 to (zeta.X_nrows-1)
){
  
  int b=a+1;
  
  int c=zeta.n_rows;
  
  while(zeta(b)==0){
    
    b=b+1;
    
    if(b==c){  //b can be higher than X.n_rows here
      break;
    }
  }
  
  return(b); //this is cpp position, not R
}


//Get the exact position below based on zeta value where zeta==1
//[[Rcpp::export]]
int lower_bound(arma::vec zeta, //zeta vector
                int a //starting point in cpp, starts from 0 to (zeta.X_nrows-1)
){
  
  
  int b=a-1;
  
  while(zeta(b)==0){
    
    if(b==0){  //b can be higher than X.n_rows here
      break;
    }
    
    b=b-1;
    
  }
  
  return(b); //this is cpp position, not R
}


//For a particular dose number get the upper bound and lower bound based on zeta value.
//NOTE:> GetBoundaries_p(1,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),c(1,0,0,0,1,1,0))
// [,1]
// [1,]  0.1
// [2,]  0.5
//But 2nd dose should be clustered with the first dose!!
//[[Rcpp::export]]
arma::vec GetBoundaries_p(int m, //initial number of dose for which we are finding boundaries
                          arma::vec p,  //toxicity vector
                          arma::vec zeta //zeta vector containing 0 and 1
){
  
  arma::vec HOLD(2);
  double bound1=0;
  double bound2=0;
  int ub=0;
  int lb=0;
  
  if(m==0){
    bound1 = 0;
    ub=upper_bound(zeta,m);
    
    if(ub==p.n_rows){
      bound2=1;
    }else{
      bound2 = p(ub);
    }
    
    HOLD(0)=bound1;
    HOLD(1)=bound2;
    
  }else{
    
    if(m==(p.n_rows-1)){
      
      lb=lower_bound(zeta,m);
      
      if(lb==0){ //1st dose
        
        if(zeta(lb)==0){  //1st dose can cluster in 0
          
          bound1=0;
          
        }else{
          
          bound1=p(lb);
        }
        
      }else{
        
        bound1=p(lb);
        
      }
      
      bound2 = 1;
      
      HOLD(0)=bound1;
      HOLD(1)=bound2;
      
    }else{  //In the middle
      
      lb=lower_bound(zeta,m);
      
      if(lb==0){
        if(zeta(lb)==0){  //1st dose can cluster in 0
          
          bound1=0;
          
        }else{
          
          bound1=p(lb);
        }
      }else{
        bound1 = p(lb);
      }
      
      ub=upper_bound(zeta,m);
      
      if(ub==zeta.n_rows){
        bound2=1;
      }else{
        bound2 = p(ub);
      }
      
      HOLD(0)=min1(bound1,bound2);
      HOLD(1)=max1(bound1,bound2);
      
    }
    
    
  }
  
  return(HOLD);
  
}


//Get a sample from truncated normal distribution whose bounds are determined by zeta vecor and p
//[[Rcpp::export]]

double TruncNormal_p(int m, //initial number of dose for which we are finding boundaries
                     arma::vec p, //toxicity vector
                     double c1,   //Std of normal distribution
                     arma::vec zeta  //zeta vector containing 0 and 1
){
  
  //m runs from 0 to (X.n_rows-1)
  //Truncation boundaries
  
  
  arma::vec BOUNDARIES = GetBoundaries_p(m,p,zeta);
  // Rprintf("bound");
  //Generate Uniform variable
  double U = as_scalar(arma::randu(1));
  
  //Get value to plug into qnorm
  // Rprintf("value");
  double X = U*R::pnorm5(BOUNDARIES(1),p(m),c1,1,0)+(1-U)*R::pnorm5(BOUNDARIES(0),p(m),c1,1,0);
  
  // Rprintf("value");
  // Rprintf("qnorm");
  X=R::qnorm5(X,p(m),c1,1,0);
  
  // Rprintf("Exit");
  
  return(X);
  
  
  
}


//[[Rcpp::export]]
double q_beta(double X, double a, double b){
  double Y=0;
  Y=R::qbeta(X,a,b,1,0);
  return(Y);
}

// [[Rcpp::export]]
double myFac(double x) { return(R::gammafn(x)); }


//[[Rcpp::export]]
double GammaFn(double z){
  //NEed some things
  double STEP = .01;
  double END=30;
  double x = 0 ; //for iterating...
  double INT_TOTAL=0;
  
  while(x<END){
    INT_TOTAL = INT_TOTAL + STEP*pow(x,z-1)*exp(-x);
    x=x+STEP;
  }
  
  return(INT_TOTAL);
}


//[[Rcpp::export]]
double BetaCDF(double z, double a, double b){
  //NEed some things
  double STEP = .001;
  double START = .001;
  double END=z;
  double x = START ; //for iterating...
  double INT_TOTAL=0;
  
  while(x<END){
    INT_TOTAL = INT_TOTAL + STEP*pow(x,a-1)*pow(1-x,b-1)*GammaFn(a+b)/(GammaFn(a)*GammaFn(b));
    x=x+STEP;
  }
  
  return(INT_TOTAL);
}


//[[Rcpp::export]]
double QBeta(double q, double a, double b){
  //The first x where F(x)>=q
  double x= .001;
  double STEP=.001;
  double PROB = 0;
  while(PROB<q){
    PROB = BetaCDF(x,a,b);
    x=x+STEP;
  }
  return(x);
  
}

//[[Rcpp::export]]

double TruncBeta(int m, //initial number of dose for which we are finding boundaries
                 arma::vec p, //toxicity vector
                 arma::vec zeta,  //zeta vector containing 0 and 1
                 double a,  //1st shape parameter of the Beta distribution
                 double b  //2nd shape parameter of the Beta distribution
){
  
  arma::vec BOUNDARIES = GetBoundaries_p(m,p,zeta);
  
  //Try writting own pbeta and qbeta functions, lower and upper bound (0.1,0.99)
  
  double Fa=R::pbeta(max1(BOUNDARIES[0],0.0001),a,b,1,0);
  double Fb=R::pbeta(min1(BOUNDARIES[1],0.9999),a,b,1,0);
  
  //,01,.99 make a difference?
  
  double U=as_scalar(arma::randu(1));
  double X=U*(Fb-Fa)+Fa;
  
  double Y=q_beta(X,a,b);
  
  return(Y);
}



//Function to randomly return a number from 1:p when it takes only p as an input

//[[Rcpp::export]]
double samp1(double p){
  
  Rcpp::IntegerVector pool = Rcpp::seq(1, p);
  
  double j=0;
  arma::vec cut_off(p);
  
  for(j=0;j<p;j++){
    cut_off(j)=(j+1)/p;
  }
  
  
  double u = as_scalar(arma::randu(1));
  
  double m=0;
  double HOLD1=m;
  for(m=0;m<p;m++){
    
    if(u<=cut_off(m)){
      
      HOLD1=m;
    }
  }
  
  return(HOLD1);
}

// [[Rcpp::export]]
Rcpp::IntegerVector VEC(double p) {
  Rcpp::IntegerVector pool = Rcpp::seq(1, p);
 // std::random_shuffle(pool.begin(), pool.end());
  return pool[Rcpp::Range(0, (p-1))];
}


////////Sample from binomial distribution/////////
//[[Rcpp::export]]
int SampBinom(int samp, double prob){
  int i;
  double U;
  int Z=0;
  for(i=0;i<samp;i++){
    
    U=as_scalar(arma::randu(1));
    
    if(U<prob){
      Z=Z+1;
    }
  }
  return(Z);
}

///////////Get the minimum dose //////////

//Revisit this!!!
//[[Rcpp::export]]
double getmin(arma::vec dose){
  double min=abs1(dose(0));
  int j=0;
  for(j=1;j<dose.n_rows;j++){
    if(dose(j)>-1000){
      if(dose(j)<min){
        min=dose(j);
      }
    }
  }
  return(min);
}

/////////Get the closest dose of MTD (either side of MTD) from a skeleton////////////


//[[Rcpp::export]]
int optdose(arma::vec dose,// skeleton
            double mtd //Maximum tolerated dose
){
  
  int m=0;
  int dosemin=0;
  arma::vec c1(dose.n_rows);
  
  //calculate the absolute value for the (dose-MTD)
  for(m=0;m<dose.n_rows;m++){
    c1(m)=abs1(dose(m)-mtd);
  }
  //get the dose which closer to MTD
  double mindose=getmin(c1);
  //Get the minimum dose level corresponding to the skeleton
  for(m=0;m<dose.n_rows;m++){
    if (c1(m)==mindose){
      dosemin=m;
    }
  }
  return(dosemin);   //add +1 to get the exact dose
  
}

///////////////////////////////////////////////////////////////////
//////////////////////////FBCRM////////////////////////////////////
///////////////////////////////////////////////////////////////////
//Function for joint posterior of p, we don't need zeta values here, as zeta_j=0 will automatically mean
//either p_j=p_(j-1) or p_(j+1) and the log contribution of X_j and Y_j will be added automatically

//[[Rcpp::export]]
double LFBCRM1( arma::vec X, //Dose assignment vector
                arma::vec Y, // Toxicity vector
                arma::vec p, //skeleton vector
                double alpha //a constant
){
  int m=0;
  
  double LogL=0;
  
  double e_alpha=pow(2.714,alpha);
  
  // Rf_PrintValue(wrap(e_alpha));
  
  for(m=0;m<X.n_rows;m++){
    
    // Rf_PrintValue(wrap(LogL));
    
    LogL = LogL+Y(m)*e_alpha*log(p(m))+(X(m)-Y(m))*log(1-pow(p(m),e_alpha));
    
    // Rf_PrintValue(wrap(Y(m)));
    // Rf_PrintValue(wrap(e_alpha*log(p(m))));
    // Rf_PrintValue(wrap(log(1-pow(p(m),e_alpha))));
    // Rf_PrintValue(wrap((X(m)-Y(m))));
    
  }
  
  return(LogL);
}


//No need to consider for doses where there are no patients allotted as those doses will cancel out in the LR when we do
//Lbeta(pi_star,zeta_star)-Lbeta(pi,zeta)
//Function for conditional posterior of p
//[[Rcpp::export]]
double Lbeta1(  arma::vec p, //skeleton vector
                arma::vec a, //vector of first parameters of beta
                arma::vec b, //vector of second parameter of beta
                arma::vec zeta
){
  
  int m=0;
  double LogL=0;
  
  for(m=0;m<p.n_rows;m++){
    if(zeta(m)==1){
      LogL = LogL+log(myFac(a(m)+b(m))/(myFac(a(m))*myFac(b(m))))+(a(m)-1)*log(p(m))+(b(m)-1)*log(1-p(m));
    }
  }
  return(LogL);
}


// [[Rcpp::export]]
double randbeta(double m, double s) { return R::rbeta(m, s); }

//[[Rcpp::export]]

double randnorm(double mu, double sigma){
  double u = as_scalar(arma::randn(1));
  double z=u*sigma+mu;
  return(z);
}

//Get the exact dose number upto which non-zero patients are allotted
//[[Rcpp::export]]
int non_zero(arma::vec X){
  int b=1;
  
  if(X(X.n_rows-1)>0){
    b=X.n_rows;
  }else{
    while(X(b)!=0){
      b=b+1;
    }
  }
  
  return(b); //this is cpp position, not R
}

//Function for conditional posterior of M (Model)

//[[Rcpp::export]]

double LM(arma::vec pi, //posterior toxicity prob. vector
          arma::vec a, //vector of first parameters of beta
          arma::vec b, //vector of second parameter of beta
          arma::vec X, //Patients allotted to each dose
          arma::vec zeta // 0/1 based on whether a dose is clustered or not
){
  
  int m=0;
  
  double LogL=0;
  
  for(m=0;m<X.n_rows;m++){
    
    // if(X(m)>0){   //Need to consider only a_beta and b_beta values for doses where X(m) is not zero
    
    if(zeta(m)==1){
      
      //Lfac doesn't work for non-integer values
      // LogL = LogL+Lfac(a(m)+b(m))-Lfac(a(m))-Lfac(b(m))+(a(m)-1)*log(pi(m))+(b(m)-1)*log(1-pi(m));
      LogL = LogL+log(myFac(a(m)+b(m))/(myFac(a(m))*myFac(b(m))))+(a(m)-1)*log(pi(m))+(b(m)-1)*log(1-pi(m));
      
    }
    
    // }
    
  }
  
  return(LogL);
}


/////////////Function for CRM likelihood/////////

//[[Rcpp::export]]
double LCRM (arma::vec X, arma::vec Y, arma::vec p, double alpha){
  
  int m=0;
  double LogL=0;
  
  for(m=0;m<Y.n_rows;m++){
    
    LogL = LogL + exp(alpha)*Y(m)*log(p(m))+(X(m)-Y(m))*log(1-pow(p(m),exp(alpha)));
  }
  
  return(exp(LogL));
}


/////////////////////////////////////////////////////////////
////////////////////////////CRM//////////////////////////////
/////////////////////////////////////////////////////////


/////////////Function for CRM likelihood for one dose/////////

//[[Rcpp::export]]
double LCRM_1D (double x, double y, double p, double al){
  
  double LogL=0;
  
  LogL =  exp(al)*y*log(p)+(x-y)*log(1-pow(p,exp(al)));
  
  return(exp(LogL));
}


////////density for normal dsitribution///////////
//[[Rcpp::export]]
double dn (double sigma, double alpha){
  
  return((1/(pow((2*3.14),0.5)*sigma))*exp(-(pow(alpha,2))/(2*pow(sigma,2))));
}


////////////Denom for updated toxicity prob calculation////////////
//[[Rcpp::export]]
double area (arma::vec X, //number of patients treated in each dose
             arma::vec Y, //Number of toxic events in each dose
             arma::vec p, //toxicity prob
             double sigma, //sd for alpha, which follows N(0,sigma^2)
             double a, //lower bound
             double b, //upper bound
             double n  //number of intervals
){
  //base of the rectangle
  double dal=(b-a)/n;
  double sum=0.0;
  //height of the rectangle
  double al=a;
  
  while (al<=b){
    sum=sum+ LCRM(X,Y,p,al) * dn(sigma,al)*dal;
    al=al+dal;
  }
  return(sum);
}


////////////Calculate the posterior probability by summing the AUC/////////////
//[[Rcpp::export]]
arma::vec areap (arma::vec X, //number of patients treated in each dose
                 arma::vec Y, //Number of toxic events in each dose
                 arma::vec p, //toxicity prob
                 double sigma, //sd for alpha, which follows N(0,sigma^2)
                 double a, //lower bound
                 double b, //upper bound
                 double n  //number of intervals
){
  //base of the rectangle
  double dal=(b-a)/n;
  int m;
  arma::vec sum(X.n_rows);
  sum.zeros();
  arma::vec Z(X.n_rows);
  Z.zeros();
  //height of the rectangle
  double c=0;
  c=area(X,Y,p,sigma,a,b,n);
  
  double al=a;
  
  while (al<=b){
    for(m=0;m<X.n_rows;m++){
      sum[m]=sum[m]+ pow(p[m],exp(al)) * LCRM(X,Y,p,al) * dn(sigma,al) * dal;
    }
    al=al+dal;
  }
  
  for(m=0;m<X.n_rows;m++){
    Z[m]=sum[m]/c;
  }
  return(Z);
}



////////////Calculate the posterior probability > MTD by summing the AUC for the first dose/////////////
//[[Rcpp::export]]
double areap_d1 (arma::vec X, //number of patients treated in each dose
                 arma::vec Y, //Number of toxic events in each dose
                 arma::vec p, //toxicity prob
                 double sigma, //sd for alpha, which follows N(0,sigma^2)
                 double mtd, //MTD toxicity for the trial
                 double a, //lower bound
                 double b, //upper bound
                 double n  //number of intervals
){
  //lower bound for the integral
  double d=0;
  d=log(log(mtd)/log(p[0]));
  //base of the rectangle
  double dal=(b-d)/n;
  double sum_1=0;
  double Z=0;
  //height of the rectangle
  double c = area(X,Y,p,sigma,a,b,n);
  b=d;
  double al=a;
  
  while (al<=b){
    
    sum_1 = sum_1 +  LCRM(X,Y,p,al) * dn(sigma,al) * dal;
    
    al=al+dal;
  }
  
  Z=sum_1/c;
  
  return(Z);
}


// Updated posterior value for p (recommended dose) after each patient cohort assignment


// [[Rcpp::export]]
List CRM_MCMC (arma::vec X, //Number of patients assigned to each dose
               arma::vec Y, // Number of toxicity events in those dose
               arma::vec p, //skeleton vector (initial toxicity level)
               double sigma, //initial variance of alpha
               double mtd,// maximum tolarated dose
               double a, //lower bound
               double b, //upper bound
               double n //number of intervals
){
  
  //Storage
  arma::vec psum(X.n_rows);
  arma::vec pbar(X.n_rows);
  
  
  //Set these equal to values
  psum.zeros();
  pbar.zeros();
  
  //posterior prob. for each dose
  psum=areap(X,Y,p,sigma,a,b,n);
  
  //optimal dose based on posterior pi values
  double opt=optdose(psum,mtd);
  
  //return the stopping prob. for 1st dose
  double L=areap_d1(X,Y,p,sigma,mtd,a,b,n);
  
  //L(D|M_k) for BMA-CRM
  double L_DMk=area(X,Y,p,sigma,a,b,n);
  
  List Z=List::create(psum,opt,L,L_DMk);
  
  return(Z);  //4 member list
}

// [[Rcpp::export]]
List CRM_RUNTRIAL(int cohort, //cohort size
                  int samp, //Maximum Sample size to run the Trial, must be divisible by cohort
                  arma::vec ptrue, //True vector of toxicity probabilities for each dose
                  arma::vec p,  //skeleton vector for toxicity
                  double sigma, //initial variance of alpha
                  double mtd,// maximum tolarated dose
                  int M, //Number of simulation
                  double p_u, //cut-off probability whether 1st dose is too toxic or not
                  double a, //lower bound
                  double b, //upper bound
                  double n  //number of intervals
){
  
  // Rprintf("Enter");
  int NCohort = samp/cohort;
  
  arma::vec X(ptrue.n_rows);
  arma::vec Y(ptrue.n_rows);
  arma::vec DoseTried(X.n_rows);
  arma::vec pbar(X.n_rows);
  
  arma::vec DoseStore(M);
  arma::mat XStore(M,ptrue.n_rows);
  arma::mat YStore(M,ptrue.n_rows);
  int i;
  int m;
  int OptDose;
  
  for(m=0;m<M;m++){
    
    OptDose=0;  //Add +1 for actual value
    X.zeros();
    Y.zeros();
    pbar.zeros();
    
    DoseTried.zeros();
    
    for(i=0;i<NCohort;i++){
      //loop for cohort
      
      X(OptDose)=X(OptDose)+cohort;
      Y(OptDose)=Y(OptDose)+SampBinom(cohort,ptrue(OptDose));
      
      DoseTried(OptDose)=1;
      
      List L=CRM_MCMC(X,Y,p,sigma,mtd,a,b,n);
      
      arma::vec pbar=L[0];
      
      double v1=L[1];
      
      double T_1=L[2];
      
      // //Create a 4 row matrix
      
      // arma::mat store(4,X.n_rows);
      // store.zeros();
      //
      // store.row(0)=X.t();
      // store.row(1)=Y.t();
      // store.row(2)=pbar.t();
      // store.row(3)=DoseTried.t();
      //
      // Rf_PrintValue(wrap(store));
      
      //best dose without any other knowledge of the trial based on only MTD
      if(sum(X)>2*cohort){
        if(T_1>p_u){
          OptDose=-2;
          break;
        }
      }
      
      if(i<(NCohort-1)){
        if(OptDose>v1){
          OptDose=OptDose-1;
        }else if(OptDose<v1){
          OptDose=OptDose+1;
        }else{
          OptDose=v1;
        }
        
      }else{
        OptDose=v1;
      }
      
    }
    
    DoseStore(m)=OptDose+1;
    XStore.row(m)=X.t();
    YStore.row(m)=Y.t();
  }
  List Z1=List::create(XStore,YStore,DoseStore);
  return(Z1);
  
}


////////////////////////////////////check on adaptive when no data
///////////////////////////////////with large data in each data
//////////////////////////////////same dataset, 5 different prosteriors with single chain
/////////////////////////////////


/////////////////////////////Calculating posterior even when dose is not tried//////////////////
//Updated posterior value for p after each patient cohort assignment


//'@importFrom Rcpp evalCpp
//'@useDynLib FBCRM
//[[Rcpp::export]]
List FBCRM_MCMC(arma::vec X, //Number of patients assigned to each dose
                arma::vec Y, // Number of toxicity events in those dose
                arma::vec p, //skeleton vector, need to calculate a_beta,b_beta
                double p_rho, //if no prior on rho then rho value, prob of not clustering
                double sigma, //Prior SD for Rho
                double mtd,// maximum tolarated dose
                int B//Number of iteration
){
  arma::vec a_beta(X.n_rows);
  arma::vec b_beta(X.n_rows);
  
  double p_theta=0.5;
  
  //Integers for loops
  int i=0;
  
  //Integers for loops
  int j=0;
  Rcpp::IntegerVector rand(X.n_rows);
  
  int m=0;
  arma::vec rho(X.n_rows);  //zeta=1 prior prob, feeding vector
  
  //dose proposal value
  arma::vec p_prop(X.n_rows);
  arma::vec p_prop1(X.n_rows);
  arma::vec p_prop2(X.n_rows);
  arma::vec zeta(X.n_rows);
  arma::vec zeta_prop(X.n_rows);
  int k=0;
  
  int l=0;
  int u=0;
  int rnum=0;
  
  arma::vec rho_prop(X.n_rows);
  
  //Storage
  arma::mat p_store(B,X.n_rows);
  arma::mat zeta_store(B,X.n_rows);
  arma::mat a_store(B,X.n_rows);
  arma::mat b_store(B,X.n_rows);
  arma::mat acc_store(B,3);
  
  arma::vec psum(X.n_rows);
  arma::vec pisum(X.n_rows);
  arma::vec zeta_sum(X.n_rows);
  arma::vec tox_tot(X.n_rows);
  arma::vec alpha_store(B);
  
  arma::vec pbar(X.n_rows);
  arma::vec pibar(X.n_rows);
  arma::vec zeta_bar(X.n_rows);
  
  pbar.zeros();
  pibar.zeros();
  zeta_bar.zeros();
  
  zeta.zeros();
  
  zeta_prop=zeta;
  
  //Set these equal to values
  psum.zeros();
  pisum.zeros();
  zeta_sum.zeros();
  tox_tot.zeros();
  p_prop.zeros();
  p_prop1.zeros();
  p_prop2.zeros();
  p_store.zeros();
  zeta_store.zeros();
  a_store.zeros();
  b_store.zeros();
  alpha_store.zeros();
  acc_store.zeros();
  
  //Used for Adaptive MH
  double B1=B;
  
  arma::vec Nprop1(X.n_rows);
  arma::vec Aprop1(X.n_rows);
  arma::vec varprop1(X.n_rows);
  
  //initialize a skeleton, a_beta and b_beta
  arma::vec a(X.n_rows);
  arma::vec b(X.n_rows);
  arma::vec a_star(X.n_rows);
  arma::vec b_star(X.n_rows);
  
  a.zeros();
  b.zeros();
  a_star.zeros();
  b_star.zeros();
  
  //starting value for skeleton matrix
  
  arma::vec RHO_STORE(B);
  arma::vec RHOVAR_STORE(B);
  
  RHO_STORE.zeros();
  RHOVAR_STORE.zeros();
  
  arma::vec z9(3);
  
  double count=X.n_rows;
  
  
  ///////////Added this, may have to remove, checking for p_rho///////////////
  int count1=non_zero(X);  //Gives the number of doses with non-zero patient allotment
  
  count1=count;
  
  //Tox counter
  double tox=0;
  
  
  double RHO=0.8;  //ESS, needed to calculate a_beta, b_beta
  double RHO_star=0;
  
  double alpha=0.01;
  double alpha_star=0.01;
  
  for(i=0;i<X.n_rows;i++){
    a_beta(i)=RHO*p(i);
    b_beta(i)=RHO*(1-p(i));
  }
  
  arma::vec a_initial=a_beta;
  arma::vec b_initial=b_beta;
  
  for(j=0;j<X.n_rows;j++){
    zeta(j)=1;    // We start with the assumption that all the doses are unclustered
    rho(j)=p_rho;  //if rho has no priors then it will take p_rho values
  }
  
  rho_prop=rho;
  
  //Used in MCMC
  double D=0;
  double U=0;
  int bound=0;
  
  //Set these equal to values
  for(m=0;m<X.n_rows;m++){
    Nprop1(m)=2;
    Aprop1(m)=1;
    varprop1(m)=.5;
  }
  
  double Da=2;
  double Dalpha=2;
  double Na=1;
  double Nalpha=2;
  
  double RHO_var=1;   //proposal distribution's variance for Rho (log-normal variance)
  double alpha_var=1;
  
  //initialize (a,b)
  for(k=0;k<X.n_rows;k++){
    a(k)=a_beta(k);
    b(k)=b_beta(k);
  }
  
  arma::vec p_use=p;
  
  arma::vec p_skel=p;
  
  
  //mixture parameter
  double theta=1;
  double theta_star=0;
  
  //moxing parameter storage vector
  arma::vec theta_store(B);
  theta_store.zeros();
  
  theta_store(0)=theta;
  
  double j1=0;
  double j2=0;
  
  for(j=0;j<B;j++){
    
    if((j%100==0) & (j<B/2)){
      
      //Adaptive variance for RHO
      
      if((Na/Da)>.6){
        //Too many accept
        RHO_var=min1((RHO_var+.1),100);
      }
      if((Na/Da)<.2){
        //Too many reject
        RHO_var=max1((RHO_var-.1),.01);
      }
      //Reset Counts
      Da=2;
      Na=1;
      
      //Adaptive proposal variance for alpha//
      
      if((Nalpha/Dalpha)>.6){
        //Too many accept
        alpha_var=alpha_var*2;
      }
      if((Nalpha/Dalpha)<.2){
        //Too many reject
        alpha_var=alpha_var/2;
      }
      //Reset Counts
      Dalpha=2;
      Nalpha=1;
      
      for(m=0;m<X.n_rows;m++){
        if((Aprop1(m)/Nprop1(m))>.5){
          //Too many accept
          varprop1(m)=min1(varprop1(m)*2,1);
        }
        if((Aprop1(m)/Nprop1(m))<.2){
          //Too many reject
          varprop1(m)=max1(varprop1(m)/2,.01);
        }
        //Reset Counts
        Nprop1(m)=2;
        Aprop1(m)=1;
      }
    }
    
    for(l=0;l<3;l++){
      z9(l)=0;
    }
    
    U=as_scalar(arma::randu(1));
    //
    // if(U<0.1){
    if(theta==1){  //it was in the FBCRM structure
      theta_star=0;  //we propose it to CRM structure
      
      //evaluate likelihood ratio, no log(alpha) or log(alpha_star) as normal distribution is
      //a symmetric distribution
      D=LFBCRM1(X,Y,p_skel,alpha) +log(p_theta)
        - LFBCRM1(X,Y,p,alpha) - log(1-p_theta)- Lbeta1(p,a,b,zeta);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        p=p_skel;
        theta=theta_star;
      }
      
    }else{  //it was CRM structure
      
      if(U<0.1){
        theta_star=1;   //propose FBCRM structure
        
        p_prop=p_use;
        
        zeta.zeros();
        zeta=zeta+1;
        
        for(j2=0;j2<10;j2++){
          for(j1=0;j1<X.n_rows;j1++){
            p_prop(j1)=TruncBeta(j1,p_prop,zeta,a_initial(j1),b_initial(j1));
            
          }
        }
        
        //a symmetric distribution
        D=LFBCRM1(X,Y,p_prop,alpha) + log(1-p_theta)+ Lbeta1(p_prop,a_initial,b_initial,zeta)
          -LFBCRM1(X,Y,p_skel,alpha) -log(p_theta);
        
        //Draw random Uniform
        U=log(as_scalar(arma::randu(1)));
        
        if(U<D){
          theta=theta_star;
          a_beta=a_initial;
          b_beta=b_initial;
          p=p_prop;
        }
      }
      
      
    }
    // }
    
    theta_store(j)=theta;
    
    if(theta==1){
      
      //Don't sample any /pi_j where X_j=0. Keep \pi_j=1
      
      //Randomly order the doses and select first m doses to be newly drawn
      rand=VEC(count);  //returns a vector
      
      for(m=1;m<count1;m++){
        rho(m)=p_rho;
      }
      
      ////////////////////alpha sampler, prior N(0,2) based on paper/////////////////////
      
      // Rprintf("Before sampling alpha");
      
      alpha_star=randnorm(alpha,alpha_var);  //proposal distribution//
      
      //evaluate likelihood ratio, no log(alpha) or log(alpha_star) as normal distribution is
      //a symmetric distribution
      D=LFBCRM1(X,Y,p,alpha_star) - LFBCRM1(X,Y,p,alpha)
        + (-pow(alpha_star,2)+pow(alpha,2))/2*pow(2,2);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        alpha=alpha_star;
        Nalpha=Nalpha+1;
      }
      
      Dalpha=Dalpha+1;
      
      alpha_store(j)=alpha;
      
      /////////////////////////
      ///joint sampler for (pi_j,zeta_j), only propose zeta=1 to zeta=0 or 0 to 1, 1st dose can't be sampled here//////////
      ////////////////////////
      
      p_prop1=p_prop2=p;
      zeta_prop=zeta;
      
      // Rprintf("Before joint sampler");
      if(count1!=1){
        
        rnum=samp1(count1-1)+1; //It will exclude first dose always     (double checked the samp function)
        
        if(rnum!=(X.n_rows-1)){  //if it's the last dose then upper bound is always 1
          bound=upper_bound(zeta,rnum);
        }else
        {
          bound=X.n_rows;
        }
        
        //1st sampler for clustering above
        if(rnum!=(X.n_rows-1)){ //if rnum is not the last dose then only clustering above can happen
          if(zeta(rnum)==1){
            if(X(rnum+1)!=0){ //Can only cluster above if the dose above has nonzero patients
              //Check if it's the last dose
              if(zeta(rnum+1)==1){  //check if /zeta_j+1=1 or not
                //Propose to cluster above
                zeta_prop(rnum+1)=0;
                p_prop1(rnum)=p_prop1(rnum+1);
                
                //Calculate joint LR but with only one change in p_prop at a time
                D = LFBCRM1(X,Y,p_prop1,alpha)+Lbeta1(p_prop1,a,b,zeta_prop)+log(1-rho(rnum))-
                  LFBCRM1(X,Y,p,alpha)- Lbeta1(p,a,b,zeta)- log(rho(rnum));
                
                //Draw random Uniform
                U=log(as_scalar(arma::randu(1)));
                
                if(U<D){
                  p=p_prop1;
                  zeta=zeta_prop;
                  z9(1)=1;
                }
              }
            }
          }
        }
        
        p_prop1=p_prop2=p;
        zeta_prop=zeta;
        
        if(rnum!=(X.n_rows-1)){  //if it's the last dose then upper bound is always 1
          bound=upper_bound(zeta,rnum);
        }else
        {
          bound=X.n_rows;
        }
        
        //2nd sampler
        if(zeta(rnum)==1){  //was unclustered, propose clustering below
          
          p_prop1(rnum)=p_prop1(rnum-1);
          
          zeta_prop(rnum)=0;
          
          //Clustering below
          if(rnum!=(X.n_rows-1)){
            //if not the last dose then replace all the values from rnum position and above with the proposal where zeta=0 until zeta=1
            for(u=rnum;u<(bound);u++){ //starting from rnum and not (rnum+1)
              p_prop1(u)=p_prop1(rnum-1);
            }
          }
          
          //Calculate joint LR but with only one change in p_prop at a time
          D = LFBCRM1(X,Y,p_prop1,alpha)+Lbeta1(p_prop1,a,b,zeta_prop)+log(1-rho(rnum))-
            LFBCRM1(X,Y,p,alpha)- Lbeta1(p,a,b,zeta)- log(rho(rnum));
          
          //Draw random Uniform
          U=log(as_scalar(arma::randu(1)));
          
          if(U<D){
            p=p_prop1;
            zeta=zeta_prop;
            z9(0)=1;
          }
          
        }
        else
        {
          //was clustered, propose unclustering
          
          zeta_prop(rnum)=1;
          
          p_use=p;
          
          p_use(rnum-1)=p_use(rnum-1)+.05;
          
          p_prop2(rnum)=TruncNormal_p(rnum,p_use,2,zeta_prop);   //May play around the variance
          
          if(rnum!=(X.n_rows-1)){
            
            for(u=(rnum+1);u<(bound);u++){
              //replace all the values above with the proposal where zeta=0 until zeta=1
              p_prop2(u)=p_prop2(rnum);
            }
          }
          
          //Calculate joint LR but with only one change in p_prop at a time, truncated normal is not symmetric
          D = LFBCRM1(X,Y,p_prop2,alpha)+ Lbeta1(p_prop2,a,b,zeta_prop)+ log(rho(rnum)) -
            LFBCRM1(X,Y,p,alpha)-Lbeta1(p,a,b,zeta)-log(1-rho(rnum));
          
          //Draw random Uniform
          U=log(as_scalar(arma::randu(1)));
          
          if(U<D){
            p=p_prop2;
            zeta=zeta_prop;
            z9(2)=1;
          }
        }
      }
      
      //Conditional sampler
      // Rprintf("Before conditional sampler");
      ///////////////////////////////////////////////////////////////////////////////////////
      //all the X.n_rows values will be proposed one at a time and a decision will be made
      //Only here 1st dose can change///////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////
      p_prop=p;
      
      for(m=0;m<count;m++){
        
        //only propose for \pi_j|\zeta_j=1
        if(zeta(rand(m)-1)==1){
          
          // p_prop(rand(m)-1) =TruncBeta((rand(m)-1),p,zeta,1+Y((rand(m)-1)),1+X(rand(m)-1)-Y(rand(m)-1));
          p_prop(rand(m)-1)=TruncNormal_p((rand(m)-1),p,varprop1((rand(m)-1)),zeta); //May play around the variance
          
          if((rand(m)-1)!=(X.n_rows-1)){ //If it's not the last dose
            bound=upper_bound(zeta,(rand(m)-1));
          }else{ //if it's the last dose then upper bound is always the last dose number
            bound=X.n_rows;
          }
          
          if((rand(m)-1)!=(X.n_rows-1)){
            for(u=rand(m);u<(bound);u++){
              //replace all the values above with the proposal where zeta=0 until zeta=1
              p_prop(u)=p_prop(rand(m)-1);
            }
          }
          //Calculate joint LR but with only one change in p_prop at a time, trun normal symmetric so no proposal ratio
          D = LFBCRM1(X,Y,p_prop,alpha)+Lbeta1(p_prop,a,b,zeta)-
            LFBCRM1(X,Y,p,alpha)- Lbeta1(p,a,b,zeta);
          
          //Draw random Uniform
          U=log(as_scalar(arma::randu(1)));
          
          
          if(U<D){
            /////////New p  accepted/////
            p=p_prop;
            Aprop1((rand(m)-1))=Aprop1((rand(m)-1))+1;
          }
          Nprop1((rand(m)-1))=Nprop1((rand(m)-1))+1;
        }
        
      }
      
      // if(j==B*.9){
      //   Rf_PrintValue(wrap(varprop1));
      //   }
      
      
      // Rf_PrintValue(wrap(p_prop));
      
      // Rprintf("Before Rho sampler");
      ///////////////////RHO (ESS) sampler///////////////////////////
      RHO_star=RHO;
      
      //randn function
      RHO_star=exp(randnorm(log(RHO),RHO_var));
      
      for(k=0;k<X.n_rows;k++){
        a_star(k)=RHO_star*p(k);
        b_star(k)=RHO_star*(1-p(k));
      }
      
      //evaluate likelihood ratio
      D=LM(p,a_star,b_star,X,zeta)+log(RHO_star) - LM(p,a,b,X,zeta)-log(RHO)
        + (-pow(RHO_star,2)+pow(RHO,2))/2*pow(sigma,2);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        RHO=RHO_star;
        a=a_star;
        b=b_star;
        Na=Na+1;
      }
      Da=Da+1;
      
      RHO_STORE[j]=RHO;
      RHOVAR_STORE[j]=RHO_var;
      
      for(k=0;k<X.n_rows;k++){
        p_store(j,k)=p(k);
      }
      
    }
    else{  //CRM structure
      
      ////////////////////alpha sampler, prior N(0,2) based on paper/////////////////////
      
      // Rprintf("Before sampling alpha");
      
      // Rf_PrintValue(wrap(p_skel));
      
      alpha_star=randnorm(alpha,alpha_var);  //proposal distribution//
      
      //evaluate likelihood ratio, no log(alpha) or log(alpha_star) as normal distribution is
      //a symmetric distribution
      D=LFBCRM1(X,Y,p_skel,alpha_star) - LFBCRM1(X,Y,p_skel,alpha);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        alpha=alpha_star;
        Nalpha=Nalpha+1;
      }
      
      Dalpha=Dalpha+1;
      
      alpha_store(j)=alpha;
      
      for(k=0;k<X.n_rows;k++){
        p_store(j,k)=p_skel(k);
      }
      
      
      // Rf_PrintValue(wrap(p_skel));
      
    }
    
    for(k=0;k<X.n_rows;k++){
      //Storage
      // p_store(j,k)=p(k);
      
      zeta_store(j,k)=zeta(k);
      
      //Storage
      a_store(j,k)=a(k);
      b_store(j,k)=b(k);
      
      //sum all the posterior p's over B/2 iterations
      if(j>=(B1/2)){
        psum(k)=psum(k)+p_store(j,k);
        pisum(k)=pisum(k)+pow(p_store(j,k),pow(2.714,alpha_store(j)));
        zeta_sum(k)=zeta_sum(k)+zeta_store(j,k);
      }
      
    }
    
    for(l=0;l<3;l++){
      acc_store(j,l)=z9(l);
    }
  }
  
  // Rprintf("After B loop");
  
  for(j=B1/2;j<B1;j++){
    
    if(pow(p_store(j,0),pow(2.714,alpha_store(j))) > mtd){
      tox=tox+1;
    }
  }
  
  for(m=0;m<X.n_rows;m++){  //make them matrix
    pbar(m)=psum(m)/(B1/2);
    pibar(m)=pisum(m)/(B1/2);
    zeta_bar(m)=zeta_sum(m)/(B1/2);
  }
  
  arma::vec mindose(1);
  
  mindose(0)=optdose(pibar,mtd);
  
  double avg_tox=tox/(B1/2);
  
  List Z=List::create(p_store,pibar,zeta_bar,mindose,zeta_store,acc_store,avg_tox,RHO_STORE,RHOVAR_STORE,a_store,b_store,rho,
                      alpha_store,pbar,theta_store);
  
  return(Z);
}



//'@importFrom Rcpp evalCpp
//'@useDynLib FBCRM
// [[Rcpp::export]]
List FBCRM_RUNTRIAL(double cohort, //cohort size
                    double max_samp, //Maximum Sample size to run the Trial, must be divisible by cohort
                    arma::vec ptrue, //True vector of toxicity probabilities for each dose
                    arma::vec p, //skeleton vector, only needed to calculate a_beta, b_beta
                    double p_rho, //if no prior on rho then rho value, zeta=1 probability
                    double sigma, //Prior SD for Rho
                    double mtd,// maximum tolarated dose
                    double p_u,// cut-off to decide if the dose is too toxic or not
                    double B,//Number of reps for each MCMC,
                    double M, //Number of simulation
                    double a, //lower bound
                    double b, //upper bound
                    double n  //number of intervals
){
  
  double NCohort = max_samp/cohort;
  
  arma::vec X(ptrue.n_rows);
  arma::vec Y(ptrue.n_rows);
  arma::vec DoseTried(X.n_rows);
  
  arma::vec DoseStore(M);
  arma::mat XStore(M,ptrue.n_rows);
  arma::mat YStore(M,ptrue.n_rows);
  arma::mat p_store(B,X.n_rows);
  int i;
  int m;
  
  // double BREAK;
  
  arma::vec OptDose(1);
  
  arma::vec pbar(ptrue.n_rows);
  arma::vec zbar(ptrue.n_rows);
  double T_1=0;
  double v1=0;
  
  pbar.zeros();
  zbar.zeros();
  
  for(m=0;m<M;m++){
    
    ///////////////////////////
    //When m=M-1, can't escalate
    ////////////////////////
    
    X.zeros();
    Y.zeros();
    OptDose.zeros();
    
    DoseTried.zeros();
    
    for(i=0;i<NCohort;i++){
      
      //loop for cohort
      
      X(OptDose(0))=X(OptDose(0))+cohort;
      
      Y(OptDose(0))=Y(OptDose(0))+SampBinom(cohort,ptrue(OptDose(0)));
      
      
      DoseTried(OptDose(0))=1;
      
      //////////////////////////////////////////////////////
      //May need to fix the p values for each trial/ Check with using the skeleton values
      ////////////////////////////////////////////////////////
      if(i<(3*NCohort/4)){
        
        List L=CRM_MCMC(X,Y,p,sigma,mtd,a,b,n);
        
        arma::vec pbar=L[0];
        
        v1=L[1];
        
        T_1=L[2];
        
      }else{
        
        List L=FBCRM_MCMC(X,Y,p,p_rho,sigma,mtd,B);
        
        v1=L[3];
        
        arma::vec pbar=L[1];
        
        //mean clustering probabilities
        arma::vec zbar=L[2];
        
        T_1=L[6];
        
      }
      
      //Create a 4 row matrix
      
      
      //best dose without any other knowledge of the trial based on only MTD
      
      //Only check the 1st dose toxicity prob iff atleast 2*cohort patients are already assigned to the 1st dose
      if(X(0)>2*cohort){
        if(T_1>p_u){
          OptDose(0)=-2;
          break;
        }
      }
      
      //If sum of dosetried is less than the total doses then we have to consider some other criteria, if not
      //optdose will be only based on FBCRM_MCMC function output
      
      if(i<(NCohort-1)){
        if(sum(DoseTried)<X.n_rows){  //if all the doses hasn't been tried
          if(OptDose(0)>v1){  //if OptDose(0)==v1, still it will increse the dose level
            OptDose=OptDose-1;
          }else{
            OptDose=OptDose+1;
          }
          
        }else{ //if all the doses has been tried
          if(OptDose(0)>v1){
            OptDose=OptDose-1;
          }else if(OptDose(0)<v1){
            OptDose=OptDose+1;
          }else{
            OptDose=v1;
          }
        }
        
        
      }else{
        OptDose=v1;
      }
      
    }
    DoseStore(m)=OptDose(0)+1;
    XStore.row(m)=X.t();
    YStore.row(m)=Y.t();
  }
  
  List Z1=List::create(XStore,YStore,DoseStore);
  return(Z1);
  
  
}


/////////////////////////////////////////////////////////////////////
////////////////////////BMA-CRM/////////////////////////////////////
////////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
List BMACRM_RUNTRIAL(int cohort, //cohort size
                     int samp, //Maximum Sample size to run the Trial, must be divisible by cohort
                     arma::vec ptrue, //True vector of toxicity probabilities for each dose
                     arma::mat skel,  //skeleton matrix for toxicity
                     double sigma, //initial variance of alpha
                     double mtd,// maximum tolarated dose
                     int M, //Number of simulation
                     double p_u,//cut off probability whether the first dose is too toxic or not
                     double a, //lower bound
                     double b, //upper bound
                     double n  //number of intervals
){
  int NCohort = samp/cohort;
  
  arma::vec X(ptrue.n_rows);
  arma::vec Y(ptrue.n_rows);
  arma::vec DoseTried(X.n_rows);
  arma::vec p(X.n_rows);
  arma::vec pbar(X.n_rows);
  arma::vec pihat(X.n_rows);
  arma::vec pstop(skel.n_rows);
  arma::vec LHOLD(skel.n_rows);
  
  arma::vec DoseStore(M);
  
  arma::mat XStore(M,ptrue.n_rows);
  arma::mat YStore(M,ptrue.n_rows);
  
  int i=0;
  int j=0;
  int k=0;
  int m=0;
  double densum;
  double toxprob;
  int OptDose=0;
  
  DoseStore.zeros();
  
  // arma::vec z9(3);
  //
  // z9.zeros();
  
  // Rf_PrintValue(wrap(XStore));
  
  
  for(m=0;m<M;m++){
    
    OptDose=0;  //Add +1 for actual value
    X.zeros();
    Y.zeros();
    DoseTried.zeros();
    p.zeros();
    pbar.zeros();
    
    for(i=0;i<NCohort;i++){
      //loop for cohort
      
      X(OptDose)=X(OptDose)+cohort;
      Y(OptDose)=Y(OptDose)+SampBinom(cohort,ptrue(OptDose));
      
      DoseTried(OptDose)=1;
      
      pihat.zeros();
      densum=0;
      toxprob=0;
      
      for(j=0;j<skel.n_rows;j++){
        
        for(k=0;k<X.n_rows;k++){
          
          p(k)=skel(j,k);
          
        }
        
        List L=CRM_MCMC(X,Y,p,sigma,mtd,a,b,n);
        
        arma::vec pbar=L[0];
        
        double T_1=L[2];
        
        pstop(j)=T_1;
        
        double L_DM=L[3];
        
        LHOLD(j)=L_DM;
        
        for(k=0;k<X.n_rows;k++){
          
          pihat(k)=pihat(k) + pbar(k) * L_DM;
          
        }
        
        densum = densum + L_DM;
        toxprob = toxprob + T_1 * L_DM;
        
      }
      
      // z9(0)=toxprob;
      // z9(1)=densum;
      
      //averaging pbars over the skeletons when all the CRMs are calculated
      pihat=pihat/densum;
      toxprob=toxprob/densum;
      
      // z9(2)=toxprob;
      //
      // Rf_PrintValue(wrap(z9));
      // Rf_PrintValue(wrap(pstop));
      // Rf_PrintValue(wrap(LHOLD));
      
      // //Create a 4 row matrix
      
      // arma::mat store(4,X.n_rows);
      // store.zeros();
      //
      // store.row(0)=X.t();
      // store.row(1)=Y.t();
      // store.row(2)=pihat.t();
      // store(3,0)=toxprob;
      //
      // Rf_PrintValue(wrap(store));
      if(sum(X)>2*cohort){
        if(toxprob>p_u){
          OptDose=-2;
          break;
        }
      }
      
      
      //If sum of dosetried is less than the total doses then we have to consider some other criteria, if not
      //optdose will be only based on CRM_MCMC function output
      
      if((sum(DoseTried))<X.n_rows){  //if all the doses has not been tried
        
        if(pihat(sum(DoseTried)-1)<mtd){
          
          OptDose=sum(DoseTried);   //in C++ this means +1
          
        }else{
          
          OptDose=optdose(pihat,mtd);
          
        }
        
      }else{  //if all the doses has been tried
        
        OptDose=optdose(pihat,mtd);
        
      }
      
      
    }
    
    DoseStore(m)=OptDose+1;
    XStore.row(m)=X.t();
    YStore.row(m)=Y.t();
    
  }
  
  // Rf_PrintValue(wrap(XStore));
  
  List Z1=List::create(XStore,YStore,DoseStore);
  return(Z1);
  
}


///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////MFBCRM////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

//Function to randomly return a number from 0:(p-1) when k is previously selected and will not be selected

//[[Rcpp::export]]

double samp2(double p, double k){
  
  double j=0;
  
  arma::vec cut_off(p-1);
  
  for(j=0;j<(p-1);j++){
    cut_off(j)=(j+1)/(p-1);
  }
  
  double u = as_scalar(arma::randu(1));
  
  double m=0;
  double HOLD1=0;
  for(m=0;m<(p-1);m++){
    
    if(u<=cut_off(m)){
      
      if(m>=k){
        HOLD1=m+1;
      }else{
        HOLD1=m;
      }
    }
  }
  
  return(HOLD1);
  
}

//Get a sample from truncated normal distribution with a upper bound of 100 and fixed variance
//[[Rcpp::export]]

double TruncNormal(double lower, //lower bound for the truncation
                   double upper, //upper bound for the truncation
                   double mean, //mean of the truncated normal
                   double var  //Variance of the truncated normal distribution
){
  
  //Generate Uniform variable
  double U = as_scalar(arma::randu(1));
  
  //Get value to plug into qnorm
  
  double X = U*R::pnorm5(upper,mean,var,1,0)+(1-U)*R::pnorm5(lower,mean,var,1,0);
  
  X=R::qnorm5(X,mean,var,1,0);
  
  return(X);
  
}



/////////////////////////////Calculating posterior even when dose is not tried//////////////////
//Updated posterior value for p after each patient cohort assignment

//'@importFrom Rcpp evalCpp
//'@useDynLib FBCRM
//[[Rcpp::export]]
List MFBCRM_MCMC(arma::vec X, //Number of patients assigned to each dose
                 arma::vec Y, // Number of toxicity events in those dose
                 arma::mat W, //skeleton matrix, need to calculate a_beta,b_beta
                 double p_rho, //if no prior on rho then rho value, prob of not clustering
                 double sigma, //Prior SD for Rho
                 double mtd,// maximum tolarated dose
                 int B//Number of iteration
){
  
  arma::mat a_beta(W.n_rows,X.n_rows);
  arma::mat b_beta(W.n_rows,X.n_rows);
  
  double p_theta=0.5;
  
  //Integers for loops
  int i=0;
  
  //Integers for loops
  int j=0;
  Rcpp::IntegerVector rand(X.n_rows);
  
  int m=0;
  arma::vec rho(X.n_rows);  //zeta=1 prior prob, feeding vector
  
  //dose proposal value
  arma::vec p_prop(X.n_rows);
  arma::vec p_prop1(X.n_rows);
  arma::vec p_prop2(X.n_rows);
  arma::vec zeta(X.n_rows);
  arma::vec zeta_prop(X.n_rows);
  int k=0;
  
  int l=0;
  int u=0;
  int rnum=0;
  
  arma::vec rho_prop(X.n_rows);
  
  //Storage
  arma::mat p_store(B,X.n_rows);
  arma::mat zeta_store(B,X.n_rows);
  arma::mat a_store(B,X.n_rows);
  arma::mat b_store(B,X.n_rows);
  arma::mat acc_store(B,3);
  
  arma::vec psum(X.n_rows);
  arma::vec pisum(X.n_rows);
  arma::vec zeta_sum(X.n_rows);
  arma::vec tox_tot(X.n_rows);
  arma::vec alpha_store(B);
  
  arma::vec pbar(X.n_rows);
  arma::vec pibar(X.n_rows);
  arma::vec zeta_bar(X.n_rows);
  
  pbar.zeros();
  pibar.zeros();
  zeta_bar.zeros();
  
  zeta.zeros();
  
  zeta_prop=zeta;
  
  //Set these equal to values
  psum.zeros();
  pisum.zeros();
  zeta_sum.zeros();
  tox_tot.zeros();
  p_prop.zeros();
  p_prop1.zeros();
  p_prop2.zeros();
  p_store.zeros();
  zeta_store.zeros();
  a_store.zeros();
  b_store.zeros();
  alpha_store.zeros();
  acc_store.zeros();
  
  //Used for Adaptive MH
  double B1=B;
  
  arma::vec Nprop1(X.n_rows);
  arma::vec Aprop1(X.n_rows);
  arma::vec varprop1(X.n_rows);
  
  //initialize a skeleton, a_beta and b_beta
  arma::vec a(X.n_rows);
  arma::vec b(X.n_rows);
  arma::vec a_star(X.n_rows);
  arma::vec b_star(X.n_rows);
  
  a.zeros();
  b.zeros();
  a_star.zeros();
  b_star.zeros();
  
  //starting value for skeleton matrix
  
  //Starting ESS value
  double RHO=0.8;  //ESS, needed to calculate a_beta, b_beta
  double RHO_star=0;
  
  arma::vec p(X.n_rows);
  
  
  //Storing the 1st row of skeleton matrix as the initial starting skeleton
  for(i=0;i<X.n_rows;i++){
    p(i)=W(0,i);
  }
  
  
  for(i=0;i<W.n_rows;i++){
    for(j=0;j<X.n_rows;j++){
      a_beta(i,j)=RHO*W(i,j);
      b_beta(i,j)=RHO*(1-W(i,j));
    }
  }
  
  arma::vec a_initial(X.n_rows);
  arma::vec b_initial(X.n_rows);
  
  arma::vec RHO_STORE(B);
  arma::vec RHOVAR_STORE(B);
  
  RHO_STORE.zeros();
  RHOVAR_STORE.zeros();
  
  
  arma::vec z9(3);
  
  double count=X.n_rows;
  
  ///////////Added this, may have to remove, checking for p_rho///////////////
  int count1=non_zero(X);  //Gives the number of doses with non-zero patient allotment
  
  count1=count;
  
  //Tox counter
  double tox=0;
  
  double alpha=0.01;
  double alpha_star=0.01;
  
  
  for(j=0;j<X.n_rows;j++){
    zeta(j)=1;    // We start with the assumption that all the doses are unclustered
    rho(j)=p_rho;  //if rho has no priors then it will take p_rho values
  }
  
  rho_prop=rho;
  
  //Used in MCMC
  double D=0;
  double U=0;
  int bound=0;
  
  //Set these equal to values
  for(m=0;m<X.n_rows;m++){
    Nprop1(m)=2;
    Aprop1(m)=1;
    varprop1(m)=.5;
  }
  
  double Da=2;
  double Dalpha=2;
  double Na=1;
  double Nalpha=2;
  
  double RHO_var=1;   //proposal distribution's variance for Rho (log-normal variance)
  double alpha_var=1;
  
  //starting value for skeleton matrix
  double HOLD=0;
  double HOLD_prop=HOLD;
  
  // Rprintf("entry 6");
  
  arma::vec HOLD_STORE(B);
  HOLD_STORE.zeros();
  
  //initialize (a,b)
  for(k=0;k<X.n_rows;k++){
    a(k)=a_beta(HOLD,k);
    b(k)=b_beta(HOLD,k);
  }
  
  arma::vec p_use=p;
  arma::vec p_skel=p;
  
  //mixture parameter
  double theta=1;
  double theta_star=0;
  
  //mixing parameter storage vector
  arma::vec theta_store(B);
  theta_store.zeros();
  
  theta_store(0)=theta;
  
  double j1=0;
  double j2=0;
  
  
  for(j=0;j<B;j++){
    
    if((j%100==0) & (j<B/2)){
      
      //Adaptive variance for RHO
      
      if((Na/Da)>.6){
        //Too many accept
        RHO_var=min1((RHO_var+.1),100);
      }
      if((Na/Da)<.2){
        //Too many reject
        RHO_var=max1((RHO_var-.1),.01);
      }
      //Reset Counts
      Da=2;
      Na=1;
      
      //Adaptive proposal variance for alpha//
      
      if((Nalpha/Dalpha)>.6){
        //Too many accept
        alpha_var=alpha_var*2;
      }
      if((Nalpha/Dalpha)<.2){
        //Too many reject
        alpha_var=alpha_var/2;
      }
      //Reset Counts
      Dalpha=2;
      Nalpha=1;
      
      for(m=0;m<X.n_rows;m++){
        if((Aprop1(m)/Nprop1(m))>.5){
          //Too many accept
          varprop1(m)=min1(varprop1(m)*2,1);
        }
        if((Aprop1(m)/Nprop1(m))<.2){
          //Too many reject
          varprop1(m)=max1(varprop1(m)/2,.01);
        }
        //Reset Counts
        Nprop1(m)=2;
        Aprop1(m)=1;
      }
    }
    
    for(l=0;l<3;l++){
      z9(l)=0;
    }
    
    U=as_scalar(arma::randu(1));
    
    if(theta==1){  //it was in the MFBCRM structure
      theta_star=0;  //we propose it to BMA-CRM structure
      
      //BMA-CRM will take the value of the initial skeleton from last iteration's chosen model
      
      for(k=0;k<X.n_rows;k++){
        p_skel(k)=W(HOLD,k);
      }
      
      //evaluate likelihood ratio, no log(alpha) or log(alpha_star) as normal distribution is
      //a symmetric distribution
      D=LFBCRM1(X,Y,p_skel,alpha) +log(p_theta)
        - LFBCRM1(X,Y,p,alpha) - log(1-p_theta)- Lbeta1(p,a,b,zeta);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        p=p_skel;
        theta=theta_star;
      }
      
    }else{  //it was BMA-CRM structure
      
      if(U<0.1){
        theta_star=1;   //propose MFBCRM structure
        
        for(k=0;k<X.n_rows;k++){
          p_prop(k)=W(HOLD,k);
          p_skel(k)=W(HOLD,k);
          a_initial(k)=a_beta(HOLD,k);
          b_initial(k)=b_beta(HOLD,k);
        }
        
        zeta.zeros();
        zeta=zeta+1;
        
        for(j2=0;j2<10;j2++){
          for(j1=0;j1<X.n_rows;j1++){
            p_prop(j1)=TruncBeta(j1,p_prop,zeta,a_initial(j1),b_initial(j1));
            
          }
        }
        
        //a symmetric distribution
        D=LFBCRM1(X,Y,p_prop,alpha) + log(1-p_theta)+ Lbeta1(p_prop,a_initial,b_initial,zeta)
          -LFBCRM1(X,Y,p_skel,alpha) -log(p_theta);
        
        //Draw random Uniform
        U=log(as_scalar(arma::randu(1)));
        
        if(U<D){
          theta=theta_star;
          for(k=0;k<X.n_rows;k++){
            a_beta(HOLD,k)=a_initial(k);
            b_beta(HOLD,k)=b_initial(k);
          }
          p=p_prop;
        }
      }
      
      
    }
    
    theta_store(j)=theta;
    
    if(theta==1){
      
      //Don't sample any /pi_j where X_j=0. Keep \pi_j=1
      
      //Randomly order the doses and select first m doses to be newly drawn
      rand=VEC(count);  //returns a vector
      
      for(m=1;m<count1;m++){
        rho(m)=p_rho;
      }
      
      ////////////////////alpha sampler, prior N(0,2) based on paper/////////////////////
      
      alpha_star=randnorm(alpha,alpha_var);  //proposal distribution//
      
      //evaluate likelihood ratio, no log(alpha) or log(alpha_star) as normal distribution is
      //a symmetric distribution
      D=LFBCRM1(X,Y,p,alpha_star) - LFBCRM1(X,Y,p,alpha)
        + (-pow(alpha_star,2)+pow(alpha,2))/2*pow(2,2);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        alpha=alpha_star;
        Nalpha=Nalpha+1;
      }
      
      Dalpha=Dalpha+1;
      
      alpha_store(j)=alpha;
      
      /////////////////////////
      ///joint sampler for (pi_j,zeta_j), only propose zeta=1 to zeta=0 or 0 to 1, 1st dose can't be sampled here//////////
      ////////////////////////
      
      p_prop1=p_prop2=p;
      zeta_prop=zeta;
      
      if(count1!=1){
        
        rnum=samp1(count1-1)+1; //It will exclude first dose always     (double checked the samp function)
        
        if(rnum!=(X.n_rows-1)){  //if it's the last dose then upper bound is always 1
          bound=upper_bound(zeta,rnum);
        }else
        {
          bound=X.n_rows;
        }
        
        //1st sampler for clustering above
        if(rnum!=(X.n_rows-1)){ //if rnum is not the last dose then only clustering above can happen
          if(zeta(rnum)==1){
            if(X(rnum+1)!=0){ //Can only cluster above if the dose above has nonzero patients
              //Check if it's the last dose
              if(zeta(rnum+1)==1){  //check if /zeta_j+1=1 or not
                //Propose to cluster above
                zeta_prop(rnum+1)=0;
                p_prop1(rnum)=p_prop1(rnum+1);
                
                //Calculate joint LR but with only one change in p_prop at a time
                D = LFBCRM1(X,Y,p_prop1,alpha)+Lbeta1(p_prop1,a,b,zeta_prop)+log(1-rho(rnum))-
                  LFBCRM1(X,Y,p,alpha)- Lbeta1(p,a,b,zeta)- log(rho(rnum));
                
                //Draw random Uniform
                U=log(as_scalar(arma::randu(1)));
                
                if(U<D){
                  p=p_prop1;
                  zeta=zeta_prop;
                  z9(1)=1;
                }
              }
            }
          }
        }
        
        p_prop1=p_prop2=p;
        zeta_prop=zeta;
        
        if(rnum!=(X.n_rows-1)){  //if it's the last dose then upper bound is always 1
          bound=upper_bound(zeta,rnum);
        }else
        {
          bound=X.n_rows;
        }
        
        //2nd sampler
        if(zeta(rnum)==1){  //was unclustered, propose clustering below
          
          p_prop1(rnum)=p_prop1(rnum-1);
          
          zeta_prop(rnum)=0;
          
          //Clustering below
          if(rnum!=(X.n_rows-1)){
            //if not the last dose then replace all the values from rnum position and above with the proposal where zeta=0 until zeta=1
            for(u=rnum;u<(bound);u++){ //starting from rnum and not (rnum+1)
              p_prop1(u)=p_prop1(rnum-1);
            }
          }
          
          //Calculate joint LR but with only one change in p_prop at a time
          D = LFBCRM1(X,Y,p_prop1,alpha)+Lbeta1(p_prop1,a,b,zeta_prop)+log(1-rho(rnum))-
            LFBCRM1(X,Y,p,alpha)- Lbeta1(p,a,b,zeta)- log(rho(rnum));
          
          //Draw random Uniform
          U=log(as_scalar(arma::randu(1)));
          
          if(U<D){
            p=p_prop1;
            zeta=zeta_prop;
            z9(0)=1;
          }
          
        }
        else
        {
          //was clustered, propose unclustering
          
          zeta_prop(rnum)=1;
          
          p_use=p;
          
          p_use(rnum-1)=p_use(rnum-1)+.05;
          
          p_prop2(rnum)=TruncNormal_p(rnum,p_use,2,zeta_prop);   //May play around the variance
          
          if(rnum!=(X.n_rows-1)){
            
            for(u=(rnum+1);u<(bound);u++){
              //replace all the values above with the proposal where zeta=0 until zeta=1
              p_prop2(u)=p_prop2(rnum);
            }
          }
          
          //Calculate joint LR but with only one change in p_prop at a time, truncated normal is not symmetric
          D = LFBCRM1(X,Y,p_prop2,alpha)+ Lbeta1(p_prop2,a,b,zeta_prop)+ log(rho(rnum)) -
            LFBCRM1(X,Y,p,alpha)-Lbeta1(p,a,b,zeta)-log(1-rho(rnum));
          
          //Draw random Uniform
          U=log(as_scalar(arma::randu(1)));
          
          if(U<D){
            p=p_prop2;
            zeta=zeta_prop;
            z9(2)=1;
          }
        }
      }
      
      //Conditional sampler
      ///////////////////////////////////////////////////////////////////////////////////////
      //all the X.n_rows values will be proposed one at a time and a decision will be made
      //Only here 1st dose can change///////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////
      p_prop=p;
      
      for(m=0;m<count;m++){
        
        //only propose for \pi_j|\zeta_j=1
        if(zeta(rand(m)-1)==1){
          
          // p_prop(rand(m)-1) =TruncBeta((rand(m)-1),p,zeta,1+Y((rand(m)-1)),1+X(rand(m)-1)-Y(rand(m)-1));
          p_prop(rand(m)-1)=TruncNormal_p((rand(m)-1),p,varprop1((rand(m)-1)),zeta); //May play around the variance
          
          if((rand(m)-1)!=(X.n_rows-1)){ //If it's not the last dose
            bound=upper_bound(zeta,(rand(m)-1));
          }else{ //if it's the last dose then upper bound is always the last dose number
            bound=X.n_rows;
          }
          
          if((rand(m)-1)!=(X.n_rows-1)){
            for(u=rand(m);u<(bound);u++){
              //replace all the values above with the proposal where zeta=0 until zeta=1
              p_prop(u)=p_prop(rand(m)-1);
            }
          }
          //Calculate joint LR but with only one change in p_prop at a time, trun normal symmetric so no proposal ratio
          D = LFBCRM1(X,Y,p_prop,alpha)+Lbeta1(p_prop,a,b,zeta)-
            LFBCRM1(X,Y,p,alpha)- Lbeta1(p,a,b,zeta);
          
          //Draw random Uniform
          U=log(as_scalar(arma::randu(1)));
          
          
          if(U<D){
            /////////New p  accepted/////
            p=p_prop;
            Aprop1((rand(m)-1))=Aprop1((rand(m)-1))+1;
          }
          Nprop1((rand(m)-1))=Nprop1((rand(m)-1))+1;
        }
        
      }
      
      ///////////////////RHO (ESS) sampler///////////////////////////
      RHO_star=RHO;
      
      //randn function
      RHO_star=exp(randnorm(log(RHO),RHO_var));
      
      for(k=0;k<X.n_rows;k++){
        a_star(k)=RHO_star*p(k);
        b_star(k)=RHO_star*(1-p(k));
      }
      
      //evaluate likelihood ratio
      D=LM(p,a_star,b_star,X,zeta)+log(RHO_star) - LM(p,a,b,X,zeta)-log(RHO)
        + (-pow(RHO_star,2)+pow(RHO,2))/2*pow(sigma,2);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        RHO=RHO_star;
        a=a_star;
        b=b_star;
        Na=Na+1;
      }
      Da=Da+1;
      
      RHO_STORE[j]=RHO;
      RHOVAR_STORE[j]=RHO_var;
      
      /////////////////////////////////////////////////////
      ///////////////Skeleton sampler//////////////////////
      ////////////////////////////////////////////////////
      HOLD_prop=HOLD;
      
      HOLD_prop=samp2(W.n_rows,HOLD);
      
      for(k=0;k<X.n_rows;k++){
        a_star(k)=RHO*W(HOLD_prop,k);
        b_star(k)=RHO*(1-W(HOLD_prop,k));
      }
      
      //Calculate joint LR with proposed model
      
      D = LM(p,a_star,b_star,X,zeta)- LM(p,a,b,X,zeta);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        
        HOLD=HOLD_prop;
        a=a_star;
        b=b_star;
      }
      
      HOLD_STORE(j)=HOLD;
      
      for(k=0;k<X.n_rows;k++){
        p_store(j,k)=p(k);
      }
      
    }
    else{  //BMACRM structure
      
      //////////////////////skeleton sampler////////////////////////////////
      
      HOLD_prop=HOLD;
      
      ///Propose to move to a new skeleton value/////////////////
      HOLD_prop=samp2(W.n_rows,HOLD);
      
      for(k=0;k<X.n_rows;k++){
        /////old skelton/////
        p(k)=W(HOLD,k);
        ////proposed skeleton/////////
        p_skel(k)=W(HOLD_prop,k);
      }
      
      //Calculate joint LR with proposed model, no a_beta and b_beta. That's why using LFBCRM1 likelihood
      
      D = LFBCRM1(X,Y,p_skel,alpha)-LFBCRM1(X,Y,p,alpha);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        HOLD=HOLD_prop;
        for(k=0;k<X.n_rows;k++){
          p(k)=p_skel(k);
        }
        
      }
      
      HOLD_STORE(j)=HOLD;
      
      ////////////////////alpha sampler, prior N(0,2) based on paper/////////////////////
      
      alpha_star=randnorm(alpha,alpha_var);  //proposal distribution//
      
      //evaluate likelihood ratio, no log(alpha) or log(alpha_star) as normal distribution is
      //a symmetric distribution
      D=LFBCRM1(X,Y,p,alpha_star) - LFBCRM1(X,Y,p,alpha);
      
      //Draw random Uniform
      U=log(as_scalar(arma::randu(1)));
      
      if(U<D){
        alpha=alpha_star;
        Nalpha=Nalpha+1;
      }
      
      Dalpha=Dalpha+1;
      
      alpha_store(j)=alpha;
      
      for(k=0;k<X.n_rows;k++){
        p_store(j,k)=p(k);
      }
      
    }
    
    for(k=0;k<X.n_rows;k++){
      
      zeta_store(j,k)=zeta(k);
      
      //Storage
      a_store(j,k)=a(k);
      b_store(j,k)=b(k);
      
      //sum all the posterior p's over B/2 iterations
      if(j>=(B1/2)){
        psum(k)=psum(k)+p_store(j,k);
        pisum(k)=pisum(k)+pow(p_store(j,k),pow(2.714,alpha_store(j)));
        zeta_sum(k)=zeta_sum(k)+zeta_store(j,k);
      }
      
    }
    
    for(l=0;l<3;l++){
      acc_store(j,l)=z9(l);
    }
  }
  
  for(j=B1/2;j<B1;j++){
    
    if(pow(p_store(j,0),pow(2.714,alpha_store(j))) > mtd){
      tox=tox+1;
    }
  }
  
  for(m=0;m<X.n_rows;m++){  //make them matrix
    pbar(m)=psum(m)/(B1/2);
    pibar(m)=pisum(m)/(B1/2);
    zeta_bar(m)=zeta_sum(m)/(B1/2);
  }
  
  arma::vec mindose(1);
  
  mindose(0)=optdose(pibar,mtd);
  
  double avg_tox=tox/(B1/2);
  
  List Z=List::create(p_store,pibar,zeta_bar,mindose,zeta_store,acc_store,avg_tox,RHO_STORE,RHOVAR_STORE,a_store,b_store,rho,
                      alpha_store,pbar,theta_store,HOLD_STORE);
  
  return(Z);
}


//'@importFrom Rcpp evalCpp
//'@useDynLib FBCRM
// [[Rcpp::export]]
List MFBCRM_RUNTRIAL(double cohort, //cohort size
                     double max_samp, //Maximum Sample size to run the Trial, must be divisible by cohort
                     arma::vec ptrue, //True vector of toxicity probabilities for each dose
                     arma::mat W, //skeleton matrix, only needed to calculate a_beta, b_beta
                     double p_rho, //if no prior on rho then rho value, zeta=1 probability
                     double sigma, //Prior SD for Rho
                     double mtd,// maximum tolarated dose
                     double p_u,// cut-off to decide if the dose is too toxic or not
                     double B,//Number of reps for each MCMC,
                     double M, //Number of simulation
                     double a, //lower bound
                     double b, //upper bound
                     double n  //number of intervals
){
  
  double NCohort = max_samp/cohort;
  
  arma::vec X(ptrue.n_rows);
  arma::vec Y(ptrue.n_rows);
  arma::vec DoseTried(X.n_rows);
  
  arma::vec DoseStore(M);
  arma::mat XStore(M,ptrue.n_rows);
  arma::mat YStore(M,ptrue.n_rows);
  arma::vec p(X.n_rows);
  arma::vec pihat(X.n_rows);
  arma::vec pstop(W.n_rows);
  arma::vec LHOLD(W.n_rows);
  arma::mat p_store(M,X.n_rows);
  int i;
  int j;
  int k;
  int m;
  double densum;
  double toxprob;
  
  arma::vec OptDose(1);
  
  arma::vec zbar(ptrue.n_rows);
  double T_1=0;
  double v1=0;
  
  zbar.zeros();
  
  for(m=0;m<M;m++){
    
    ///////////////////////////
    //When m=M-1, can't escalate
    ////////////////////////
    
    X.zeros();
    Y.zeros();
    
    OptDose.zeros();
    
    DoseTried.zeros();
    
    
    for(i=0;i<NCohort;i++){
      
      //loop for cohort
      
      X(OptDose(0))=X(OptDose(0))+cohort;
      
      Y(OptDose(0))=Y(OptDose(0))+SampBinom(cohort,ptrue(OptDose(0)));
      
      
      DoseTried(OptDose(0))=1;
      
      
      //////////////////////////////////////////////////////
      //May need to fix the p values for each trial/ Check with using the skeleton values
      ////////////////////////////////////////////////////////
      if(i<(3*NCohort/4)){
        
        pihat.zeros();
        densum=0;
        toxprob=0;
        
        for(j=0;j<W.n_rows;j++){
          
          for(k=0;k<X.n_rows;k++){
            
            p(k)=W(j,k);
            
          }
          
          List L=CRM_MCMC(X,Y,p,sigma,mtd,a,b,n);
          
          arma::vec pbar=L[0];
          
          T_1=L[2];
          
          pstop(j)=T_1;
          
          double L_DM=L[3];
          
          LHOLD(j)=L_DM;
          
          for(k=0;k<X.n_rows;k++){
            
            pihat(k)=pihat(k) + pbar(k) * L_DM;
            
          }
          
          densum = densum + L_DM;
          toxprob = toxprob + T_1 * L_DM;
          
        }
        
        //averaging pbars over the skeletons when all the CRMs are calculated
        pihat=pihat/densum;
        
        toxprob=toxprob/densum;
        
        T_1=toxprob;
        
        v1=optdose(pihat,mtd);
        
      }else{
        
        List L=MFBCRM_MCMC(X,Y,W,p_rho,sigma,mtd,B);
        
        //Optimal dose
        v1=L[3];
        
        //mean clustering probabilities
        arma::vec zbar=L[2];
        
        //Stopping probability for the 1st dose
        T_1=L[6];
        
        
      }
      
      
      //best dose without any other knowledge of the trial based on only MTD
      
      //Only check the 1st dose toxicity prob iff atleast 2*cohort patients are already assigned to the 1st dose
      if(X(0)>2*cohort){
        if(T_1>p_u){
          OptDose(0)=-2;
          break;
        }
      }
      
      //If sum of dosetried is less than the total doses then we have to consider some other criteria, if not
      //optdose will be only based on FBCRM_MCMC function output
      
      if(i<(NCohort-1)){
        if(sum(DoseTried)<X.n_rows){  //if all the doses hasn't been tried
          if(OptDose(0)>v1){  //if OptDose(0)==v1, still it will increse the dose level
            OptDose=OptDose-1;
          }else{
            OptDose=OptDose+1;
          }
          
        }else{ //if all the doses has been tried
          if(OptDose(0)>v1){
            OptDose=OptDose-1;
          }else if(OptDose(0)<v1){
            OptDose=OptDose+1;
          }else{
            OptDose=v1;
          }
        }
        
        
      }else{
        OptDose=v1;
      }
      
      
    }
    
    DoseStore(m)=OptDose(0)+1;
    XStore.row(m)=X.t();
    YStore.row(m)=Y.t();
    
    
  }
  
  List Z1=List::create(XStore,YStore,DoseStore);
  return(Z1);
  
  
}




