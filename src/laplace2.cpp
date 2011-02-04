//#include "models.h"
#include "utils.h"


RcppExport SEXP nsem2(SEXP data,
		      SEXP theta,
		      SEXP Sigma,
		      SEXP modelpar,
		      SEXP control);

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
    h(pos) = data(i)-mu1(i)-lambda1(i)*eta(0);
  }
  for (unsigned i=0; i<mu2.n_elem; i++) {
    pos++;
    h(pos) = data(i+mu1.n_elem)-mu2(i)-lambda2(i)*eta(1);
  }
  h(1) = eta(1)-gamma(0)*eta(0)-gamma(1)*eta(0)*eta(0);
  h(0) = eta(0);
  //  cerr << "her\n";
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
  //  D[1,1] = 1; D[2,2] = 1;
  for (unsigned i=0; i<mu1.n_elem; i++) {
    D(i+2,0) = -lambda1(i);
}
  for (unsigned i=0; i<mu2.n_elem; i++) {
    D(i+2+mu1.n_elem,1) = -lambda2(i);
  }
  // cerr << "D=\n" << D << endl;
  mat dS = -trans(D)*iS;  
  res.grad = dS*h;
  res.hess = dS*D;
  //mat dummy = trans(iSh)*d2eta0;
  //res.hess(0,0) -= dummy[0];
  res.hess(1,0) -= -2*gamma(1)*iSh(2);
  return(res);
}

rowvec laNR2(const rowvec &data, const mat &iS, const double &detS,
	     const rowvec &mu1, const rowvec &mu2,
	     const rowvec &lambda1, const rowvec &lambda2,
	     const rowvec &beta1,
	     const rowvec &beta2,
	     const rowvec &gamma,
	     const double &Dtol, const unsigned &niter, const double &lambda) {  
  rowvec eta = zeros(1,2);
  for (unsigned i=0; i<niter; i++) {
    hObj K = h2(eta,data,iS,
		mu1,mu2,lambda1,lambda2,beta1,beta2,gamma);
      double Sabs = as_scalar(trans(K.grad)*K.grad);      
      if (Sabs<Dtol) break;
      //      mat Delta = trans(K.grad)*inv(K.hess + 0.1*eye(K.hess.n_cols,K.hess.n_cols));
      mat Delta = trans(K.grad)*inv(K.hess);
      eta = eta-lambda*Delta;  
      //      hObj K = h(eta1,data,iS,
      //		mu1,mu2,lambda1,lambda2,beta,gamma);
      
  }

  hObj K = h2(eta,data,iS,
	      mu1,mu2,lambda1,lambda2,beta1,beta2,gamma);
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
  rowvec res(3);
  res(0) = logI;  for (unsigned i=0; i<2; i++) res(i+1) = eta(i);
  return(res);
}

RcppExport SEXP nsem2(SEXP data,  
		      SEXP theta,
		      SEXP Sigma,
		      SEXP modelpar,
		      SEXP control
		      ) {   

  srand ( time(NULL) ); /* initialize random seed: */
  
  Rcpp::NumericVector Theta(theta);  
  Rcpp::NumericMatrix D(data);
  unsigned nobs = D.nrow(), k = D.ncol();
  mat Data(D.begin(), nobs, k, false); // Avoid copying
  Rcpp::NumericMatrix V(Sigma);  
  mat S(V.begin(), V.nrow(), V.ncol()); 
  mat iS = inv(S);
  double detS = det(S);

  RcppParams Modelpar(modelpar);
  unsigned nlatent = Modelpar.getIntValue("nlatent");
  unsigned ny1 = Modelpar.getIntValue("nvar1");
  unsigned ny2 = Modelpar.getIntValue("nvar2");
  unsigned npred1 = Modelpar.getIntValue("npred1");
  unsigned npred2 = Modelpar.getIntValue("npred2");
  RcppParams Control(control);   
  double lambda = Control.getDoubleValue("lambda");
  double niter = Control.getIntValue("niter");
  double Dtol = Control.getDoubleValue("Dtol");

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

  // rowvec eta(nlatent);
  // eta[0]=1; eta[1]=2; 
  // cerr << "mu1=" << mu1 << endl;
  // cerr << "mu2=" << mu2 << endl;
  // cerr << "l1=" << lambda1 << endl;
  // cerr << "l2=" << lambda2 << endl;
  // cerr << "beta1=" << beta1 << endl;
  // cerr << "beta2=" << beta2 << endl;
  // cerr << "gamma=" << gamma << endl;  
  // cerr << "invSigma=\n" << iS << endl;
  
  mat lap(nobs,3);
  for (unsigned i=0; i<nobs; i++) {
    rowvec newlap = laNR2(Data.row(i), iS, detS,
			  mu1, mu2, lambda1, lambda2, beta1, beta2, gamma,
			 Dtol,niter,lambda);
    lap.row(i) = newlap;
  }

  List  res;
  int dimeta = 2;
  res["indiv"] = lap;
  res["logLik"] = sum(lap.col(0)) + (dimeta-V.nrow())*log(2.0*math::pi())*nobs/2;
  res["norm0"] = (dimeta-V.nrow())*log(2*math::pi())/2;
  return res;
}



