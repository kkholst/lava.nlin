//#include "models.h"
#include "utils.h"
#include "quadrule.h"
#include <cmath>

/*
c(y1,y2,y3,...) ~ eta2
c(x1,x2,x3,...) ~ eta1
eta2 ~ eta1 + eta1^2
c(eta2,eta1) ~ x
*/
RcppExport SEXP nsem2(SEXP data,
		      SEXP theta,
		      SEXP Sigma,
		      SEXP modelpar,
		      SEXP control,
		      SEXP eb0);

hObj h2(const rowvec &eta,
	const rowvec &data, const mat &iS,
	const rowvec &mu1, const rowvec &mu2,
	const rowvec &lambda1, const rowvec &lambda2,
	const rowvec &beta1,
	const rowvec &beta2,
	const rowvec &gamma,
	int fast=0) {  

  unsigned ny=mu1.n_elem+mu2.n_elem;
  int k=2+ny;
  colvec h(k);
  unsigned pos=1;
  for (unsigned i=0; i<mu1.n_elem; i++) {
    pos++;    
    h(pos) = data(i)-mu1(i)*(i>0)-lambda1(i)*eta(0);
  }
  for (unsigned i=0; i<mu2.n_elem; i++) {
    pos++;
    h(pos) = data(i+mu1.n_elem)-mu2(i)*(i>0)-lambda2(i)*eta(1);
  }
  h(1) = eta(1)-gamma(0)*eta(0)-gamma(1)*eta(0)*eta(0);
  h(0) = eta(0);
  h(1) -= mu2(0);
  h(0) -= mu1(0);
  unsigned dpos = ny-1;
  for (unsigned i=0; i<beta1.n_elem; i++) {
    dpos++;
    h(0) -= beta1(i)*data(dpos);
  }
  for (unsigned i=0; i<beta2.n_elem; i++) {
    dpos++;
    h(1) -= beta2(i)*data(dpos);
  }
  mat iSh = iS*h;
  hObj res;
  res.h = h;
  res.hSh = -0.5*as_scalar(trans(h)*iSh);
  if (fast==1) 
    return(res);

  mat D = zeros(k,2);
  D(0,0) = 1;  D(1,0) = -gamma(0)-2*gamma(1)*eta(0);
  D(1,1) = 1; 
  for (unsigned i=0; i<mu1.n_elem; i++) {
    D(i+2,0) = -lambda1(i);
}
  for (unsigned i=0; i<mu2.n_elem; i++) {
    D(i+2+mu1.n_elem,1) = -lambda2(i);
  }
  mat dS = -trans(D)*iS;  
  res.grad = dS*h;
  res.hess = dS*D;
  res.hess(1,0) -= -2*gamma(1)*iSh(2);
  return(res);
}

rowvec laNR2(rowvec eta,
	     const rowvec &data, const mat &iS, const double &detS,
	     const rowvec &mu1, const rowvec &mu2,
	     const rowvec &lambda1, const rowvec &lambda2,
	     const rowvec &beta1,
	     const rowvec &beta2,
	     const rowvec &gamma,
	     const double &Dtol, const unsigned &niter, const double &lambda, const int &nq=0) {


  double ll=0;
  for (unsigned i=0; i<niter; i++) {
    hObj K = h2(eta,data,iS,
		    mu1,mu2,lambda1,lambda2,beta1,beta2,gamma);
    ll = K.hSh;
    double Sabs = norm(K.grad,2);
    if (Sabs<Dtol) break;
    mat Delta = trans(K.grad)*inv(K.hess);
    eta = eta-lambda*Delta;
  }

  rowvec eta0(1,2);
  eta0(0) = mu1(0);
  eta0(1) = mu2(0);
  /*
  double diff = norm(eta0-eta,2);
  if (diff>0.5) { // Try starting the optimization in the mean eta0
    double ll0=0;
      for (unsigned i=0; i<niter; i++) {
	hObj K = h2(eta0,data,iS,
		    mu1,mu2,lambda1,lambda2,beta1,beta2,gamma);
	ll0 = K.hSh;
	double Sabs = norm(K.grad,2);
	if (Sabs<Dtol) break;
	mat Delta = trans(K.grad)*inv(K.hess);
	eta0 = eta0-lambda*Delta;      
      }
      if (ll0>ll) eta=eta0;
  }
*/
    
  unsigned p = (iS.n_cols-2);
  rowvec res(3);
  hObj K = h2(eta,data,iS,
	      mu1,mu2,lambda1,lambda2,beta1,beta2,gamma);
  double logHdet;
  double sign;
  if (nq==0) {
    log_det(logHdet,sign,K.hess); // log(det(-K.hess))
    if (std::isnan(logHdet)) logHdet = -1000;
    double logI = K.hSh - 0.5*(logHdet+log(detS));
    res(0) = logI -p*0.5*log(2*datum::pi);
  } else {

    QuadRule gh(nq);
    vec z=gh.Abscissa();
    vec w=gh.Weight();
    double C = sqrt(2);
    mat H = K.hess;

    bool useSVD=false;
    mat G = -H;
    double logdetG=0;
    mat Bi;
    if (!useSVD)  {
      mat U; vec s; mat V;
      svd(U,s,V,G);
      vec is=s;
      for (unsigned k=0; k<s.n_elem; k++) {
	if (s[k]<1e-9) {
	  s[k] = 1e-9;
	  is[k] = 0;
	} else {
	  is[k] = 1/sqrt(s[k]);
	}
      }
      Bi = trans(U*diagmat(is)*V);
      logdetG = sum(log(s));
    } else {
      mat B = chol(G);
      Bi = B.i();
      logdetG = log(det(G));
    }

    double Sum = 0;
    for (unsigned k=0; k<z.n_elem; k++) {
      for (unsigned l=0; l<z.n_elem; l++) {
	mat z0(2,1);
	z0(0) = z[k]; z0(1) = z[l];
	rowvec a0 = eta+C*trans(Bi*z0);
	//rowvec a0 = eta+C*trans(B.i()*z0);
	K = h2(a0,data,iS,mu1,mu2,lambda1,lambda2,beta1,beta2,gamma,true);
	double w0 = w[k]*w[l]*exp(z0[0]*z0[0])*exp(z0[1]*z0[1]);
 	double ll0 = -0.5*log(detS) + K.hSh - (p+2)*0.5*log(2*datum::pi);
	Sum += exp(ll0)*w0;
      }
    }
    res(0) = 2*log(C)-0.5*logdetG+log(Sum);
    // res(0) = 2*log(C)-0.5*log(det(G))+log(Sum);    
      
  }
  for (unsigned i=0; i<2; i++) res(i+1) = eta(i);
  return(res); 
}

RcppExport SEXP nsem2(SEXP data,  
		      SEXP theta,
		      SEXP Sigma,
		      SEXP modelpar,
		      SEXP control,
		      SEXP eb0
		      ) {   
BEGIN_RCPP
  //  srand ( time(NULL) ); /* initialize random seed: */
  
  Rcpp::NumericVector Theta(theta);  
  Rcpp::NumericMatrix D(data);
  unsigned nobs = D.nrow(), k = D.ncol();
  mat Data(D.begin(), nobs, k, false); // Avoid copying
  Rcpp::NumericMatrix V(Sigma);  
  mat S(V.begin(), V.nrow(), V.ncol()); 
  mat iS = inv(S);
  double detS = det(S);  
  Rcpp::NumericMatrix EB0(eb0);
  mat EBstart(EB0.begin(),EB0.nrow(),EB0.ncol());

  Rcpp::List Modelpar(modelpar);
  Rcpp::IntegerVector _ny1 = Modelpar["nvar1"]; unsigned ny1 = _ny1[0];
  Rcpp::IntegerVector _ny2 = Modelpar["nvar2"]; unsigned ny2 = _ny2[0];
  Rcpp::IntegerVector _npred1 = Modelpar["npred1"]; unsigned npred1 = _npred1[0];
  Rcpp::IntegerVector _npred2 = Modelpar["npred2"]; unsigned npred2 = _npred2[0];
  Rcpp::List Control(control);   
  Rcpp::NumericVector _lambda = Control["lambda"]; double lambda = _lambda[0];
  Rcpp::NumericVector _niter = Control["niter"]; double niter = _niter[0];
  Rcpp::NumericVector _Dtol = Control["Dtol"]; double Dtol = _Dtol[0];
  Rcpp::IntegerVector _nq = Control["nq"]; unsigned nq = _nq[0];

  rowvec mu1(ny1), lambda1(ny1);
  rowvec mu2(ny2), lambda2(ny2);
  rowvec beta1(npred1); 
  rowvec beta2(npred2); 
  rowvec gamma(2);
  unsigned pos=0;
  for (unsigned i=0; i<ny1; i++) {
    mu1(i) = Theta[pos];
    pos++;
  }
  for (unsigned i=0; i<ny2; i++) {
    mu2(i) = Theta[pos];
    pos++;
  }
  lambda1(0) = 1;
  for (unsigned i=1; i<ny1; i++) {
    lambda1(i) = Theta[pos];
    pos++;
  }
  lambda2(0) = 1;
  for (unsigned i=1; i<ny2; i++) {
    lambda2(i) = Theta[pos];
    pos++;
  }
  for (unsigned i=0; i<npred1; i++) {
    beta1(i) = Theta[pos];
    pos++;
  }
  for (unsigned i=0; i<npred2; i++) {
    beta2(i) = Theta[pos];
    pos++;
  }
  gamma(0) = Theta[pos]; gamma(1) = Theta[pos+1];

  rowvec eta(1,2); eta(0)=mu1(0); eta(1)=mu2(0);

  mat lap(nobs,3);
  for (unsigned i=0; i<nobs; i++) {
    if (EBstart.n_rows==nobs) eta = EBstart.row(i);
    // Rcpp::Rcout << "eta=" << eta << endl;
    rowvec newlap = laNR2(eta,
			  Data.row(i), iS, detS,
			  mu1, mu2, lambda1, lambda2, beta1, beta2, gamma,
			  Dtol,niter,lambda,nq);
    lap.row(i) = newlap;
  }

  List  res;
  res["indiv"] = lap;
  res["logLik"] = sum(lap.col(0));
  //    (dimeta-V.nrow())*log(2.0*datum::pi)*nobs/2;
  // res["norm0"] = (dimeta-V.nrow())*log(2*datum::pi)/2;
  return res;
END_RCPP  
}



