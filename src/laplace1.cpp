//#include "models.h"
#include "utils.h"
//#include <time.h>

RcppExport SEXP nsem3(SEXP data,
				SEXP theta,
				SEXP Sigma,
				SEXP modelpar,
				SEXP control);
				
hObj h(const rowvec &eta,
       const rowvec &data, const mat &iS,
       const rowvec &mu0, const rowvec &mu1, const rowvec &mu2,
       const rowvec &lambda0, 
       const rowvec &lambda1, 
       const rowvec &lambda2,
       const rowvec &beta0,
       const rowvec &beta1,
       const rowvec &beta2,
       const rowvec &gamma,
       int fast=0) {  

  unsigned ny=mu0.n_elem+mu1.n_elem+mu2.n_elem;
  int k=3+ny;
  colvec h(k);
  unsigned pos=2;
  unsigned dpos = 0;
  //  cerr << "mu0" << endl;
  for (unsigned i=0; i<mu0.n_elem; i++) {
    pos++;
    //    cerr << data(dpos) << "; ";
    h(pos) = data(dpos)-mu0(i)-lambda0(i)*eta(0);
    dpos++;
  }
  //  cerr << "mu1" << endl;
  for (unsigned i=0; i<mu1.n_elem; i++) {
    pos++;
    //    cerr << data(dpos) << "; ";
    h(pos) = data(dpos)-mu1(i)-lambda1(i)*eta(1);
    dpos++;
  }
  //  cerr << "mu2" << endl;
  for (unsigned i=0; i<mu2.n_elem; i++) {
    pos++;
    //    cerr << data(dpos) << "; ";
    h(pos) = data(dpos)-mu2(i)-lambda2(i)*eta(2);
    dpos++;
  }
  double zeta2 = eta(2)-gamma(0)*eta(0)-gamma(1)*eta(0)*eta(0);
  double zeta1 = eta(1)-eta(0);
  double zeta0 = eta(0);
  // cerr << "dpos=" << dpos << endl;
  // cerr << "mu0=" << mu0 << endl;
  // cerr << "mu1=" << mu1 << endl;
  // cerr << "mu2=" << mu2 << endl;
  //  cerr << "beta" << endl;
  for (unsigned i=0; i<beta0.n_elem; i++) {
    //    cerr << data(dpos) << "; ";
    zeta0 -= beta0(i)*data(dpos);
    dpos++;
  }
  for (unsigned i=0; i<beta1.n_elem; i++) {
    zeta1 -= beta1(i)*data(dpos);
    dpos++;
  }
  for (unsigned i=0; i<beta2.n_elem; i++) {
    zeta2 -= beta2(i)*data(dpos);
    dpos++;
  }
  h(0) = zeta0;
  h(1) = zeta1;
  h(2) = zeta2;

  mat iSh = iS*h;
  hObj res;
  res.h = h;
  res.hSh = -0.5*as_scalar(trans(h)*iSh);
  if (fast==1) 
    return(res);
  
  mat D = zeros(k,3);
  D(0,0) = 1;  D(1,0) = -1; D(2,0) = -gamma(0)-2*gamma(1)*eta(0);
  D(1,1) = 1; D(2,2) = 1;
  //  D[1,1] = 1; D[2,2] = 1;
  pos = 3;
  for (unsigned i=0; i<mu0.n_elem; i++) {
    D(pos,0) = -lambda0(i);
    pos++;
  }
  for (unsigned i=0; i<mu1.n_elem; i++) {
    D(pos,1) = -lambda1(i);
    pos++;
  }
  for (unsigned i=0; i<mu2.n_elem; i++) {
    D(pos,2) = -lambda2(i);
    pos++;
  }
  // cerr << "D=\n" << D << endl;
  mat dS = -trans(D)*iS;  
  res.grad = dS*h;
  res.hess = dS*D;
  //mat dummy = trans(iSh)*d2eta0;
  //res.hess(0,0) -= dummy[0];
  res.hess(0,0) -= -2*gamma(1)*iSh(2);
  return(res);
}

rowvec laNR(const rowvec &data, const mat &iS, const double &detS,
	    const rowvec &mu0, const rowvec &mu1, const rowvec &mu2,
	    const rowvec &lambda0, const rowvec &lambda1, const rowvec &lambda2,
	    const rowvec &beta0, const rowvec &beta1, const rowvec &beta2,
	    const rowvec &gamma,
	    const double &Dtol, const unsigned &niter, const double &lambda) {  
  rowvec eta = zeros(1,3);
  for (unsigned i=0; i<niter; i++) {
      hObj K = h(eta,data,iS,
		 mu0,mu1,mu2,lambda0,lambda1,lambda2,beta0,beta1,beta2,gamma);
      double Sabs = as_scalar(trans(K.grad)*K.grad);      
      if (Sabs<Dtol) break;
      //      mat Delta = trans(K.grad)*inv(K.hess + 0.1*eye(K.hess.n_cols,K.hess.n_cols));
      mat Delta = trans(K.grad)*inv(K.hess);
      eta = eta-lambda*Delta;  
      //      hObj K = h(eta1,data,iS,
      //		mu1,mu2,lambda1,lambda2,beta,gamma);
      
  }

  hObj K = h(eta,data,iS,
	     mu0,mu1,mu2,lambda0,lambda1,lambda2,beta0,beta1,beta2,gamma);
  //  cerr << "K.grad=" << K.grad << endl;
  double logHdet;
  double sign;
  log_det(logHdet,sign,K.hess); // log(det(-K.hess))
  if (isnan(logHdet)) logHdet = -1000;
  double logI = K.hSh - 0.5*(logHdet+log(detS));
  //  double logI = K.hSh - 0.5*log(detS);
  //  cerr << "logI" << logI << endl;
  //  cerr << "hess" << -K.hess << endl;
  //  cerr << "hSh" << -K.hSh << endl;
  rowvec res(4);
  res(0) = logI;  for (unsigned i=0; i<3; i++) res(i+1) = eta(i);
  return(res);
}

RcppExport SEXP nsem3(SEXP data,  
		      SEXP theta,
		      SEXP Sigma,
		      SEXP modelpar,
		    SEXP control
		    ) {   

  //  srand ( time(NULL) ); /* initialize random seed: */
  
  Rcpp::NumericVector Theta(theta);  
  Rcpp::NumericMatrix D(data);
  unsigned nobs = D.nrow(), k = D.ncol();
  mat Data(D.begin(), nobs, k, false); // Avoid copying
  Rcpp::NumericMatrix V(Sigma);  
  mat S(V.begin(), V.nrow(), V.ncol()); 
  mat iS = inv(S);
  double detS = det(S);
 


  Rcpp::List Modelpar(modelpar);
  Rcpp::IntegerVector _nlatent = Modelpar["nlatent"]; unsigned nlatent = _nlatent[0];
  Rcpp::IntegerVector _ny0 = Modelpar["nvar0"]; unsigned ny0 = _ny0[0];
  Rcpp::IntegerVector _ny1 = Modelpar["nvar1"]; unsigned ny1 = _ny1[0];
  Rcpp::IntegerVector _ny2 = Modelpar["nvar2"]; unsigned ny2 = _ny2[0];
  Rcpp::IntegerVector _npred0 = Modelpar["npred0"]; unsigned npred0 = _npred0[0];
  Rcpp::IntegerVector _npred1 = Modelpar["npred1"]; unsigned npred1 = _npred1[0];
  Rcpp::IntegerVector _npred2 = Modelpar["npred2"]; unsigned npred2 = _npred2[0];
  Rcpp::List Control(control);   
  Rcpp::NumericVector _lambda = Control["lambda"]; double lambda = _lambda[0];
  Rcpp::NumericVector _niter = Control["niter"]; double niter = _niter[0];
  Rcpp::NumericVector _Dtol = Control["Dtol"]; double Dtol = _Dtol[0];


  rowvec mu0(ny0), lambda0(ny0);
  rowvec mu1(ny1), lambda1(ny1);
  rowvec mu2(ny2), lambda2(ny2);
  rowvec beta0(npred0); 
  rowvec beta1(npred1); 
  rowvec beta2(npred2);
  rowvec gamma(2);
  unsigned pos=0;
  for (unsigned i=0; i<ny0; i++) {
    mu0(i) = Theta[pos];
    pos++;
  }
  for (unsigned i=0; i<ny1; i++) {
    mu1(i) = Theta[pos];
    pos++;
  }
  for (unsigned i=0; i<ny2; i++) {
    mu2(i) = Theta[pos];
    pos++;
  }
  for (unsigned i=0; i<ny0; i++) {
    lambda0(i) = Theta[pos];
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
  for (unsigned i=0; i<npred0; i++) {
    beta0(i) = Theta[pos];
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

  // cerr << "mu0=" << mu0 << endl;
  // cerr << "mu1=" << mu1 << endl;
  // cerr << "mu2=" << mu2 << endl;
  // cerr << "lambda0=" << lambda0 << endl;
  // cerr << "lambda1=" << lambda1 << endl;
  // cerr << "lambda2=" << lambda2 << endl;
  // cerr << "beta0=" << beta0 << endl;
  // cerr << "beta1=" << beta1 << endl;
  // cerr << "beta2=" << beta2 << endl;
  // cerr << "gamma=" << gamma << endl;
  
  mat lap(nobs,4);
  for (unsigned i=0; i<nobs; i++) {
    rowvec newlap = laNR(Data.row(i), iS, detS,
			 mu0, mu1, mu2, 
			 lambda0, lambda1, lambda2, 
			 beta0,beta1, beta2, gamma,
			 Dtol,niter,lambda);
    lap.row(i) = newlap;
  }

  List  res;
  res["indiv"] = lap;
  res["logLik"] = sum(lap.col(0)) + (3-V.nrow())*log(2.0*math::pi())*nobs/2;
  res["norm0"] = (3-V.nrow())*log(2*math::pi())/2;
  return res;
}



