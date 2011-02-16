#ifndef MODELS_H
#define MODELS_H

/* ************************************************ */

/* {{{ Typedefs and prototypes...*/

typedef double (*DesignFunPt)(string model, const colvec &theta,
			      const mat &y,const mat &eta, int clustersize,
			      const SEXP &modelpar);
typedef double (DesignFun)(string model, const colvec &theta,
			   const mat &y,const mat &eta, int clustersize,
			   const SEXP &modelpar);

int lookup(const char *model);
DesignFunPt DensFixed(string model);

/* }}} */

/* ************************************************ */

/* {{{ Model prototypes & DensFixed function...*/

DesignFun evalh;
DesignFun weibull;
DesignFun weibull2;
DesignFun weibullmm;
DesignFun weibullmm2;
DesignFun weibullmm3;
DesignFun mimic;
DesignFun mimic2;
DesignFun nsem1;

const char* lookup_table[] = { 			       
  "weibull",
  "weibull2",
  "weibullmm",
  "weibullmm2",
  "weibullmm3",
  "mimic",
  "mimic2",
  "nsem1"
};
enum { ITEM_WEIBULL,
       ITEM_WEIBULL2,
       ITEM_WEIBULLMM,
       ITEM_WEIBULLMM2,
       ITEM_WEIBULLMM3,
       ITEM_MIMIC,
       ITEM_MIMIC2,
       ITEM_NSEM1,
       ITEM_NONE
};

DesignFunPt DensFixed(string model) {
  const char *mymodel = model.c_str();
  switch (lookup(mymodel)) {
  case ITEM_WEIBULL:
    return(weibull);
    break;
  case ITEM_WEIBULL2:
    return(weibull2);
    break;
  case ITEM_WEIBULLMM:
    return(weibullmm);
    break;
  case ITEM_WEIBULLMM2:
    return(weibullmm2);
    break;
  case ITEM_WEIBULLMM3:
    return(weibullmm3);
    break;
  case ITEM_MIMIC:
    return(mimic);
    break;
  case ITEM_MIMIC2:
    return(mimic2);
    break;
  case ITEM_NSEM1:
      return(nsem1);
      break;
  default : 
    return(evalh);
  }
  return(evalh);
}

/* case ITEM_PROBIT1: { */
/*   Matrix<double, Col> NewS(1,1); */
/*   return(PROBIT1); */
/*   break; */
/* } */

/* }}} */

/* ************************************************ */

/* {{{ evalh, lookup...*/

int lookup(const char model[]) {
  const int available_models = sizeof lookup_table / sizeof *lookup_table;
  for(int i=0; i!= available_models; i++)
    if( strcmp( model, lookup_table[i]) == 0 )
      return i;
  return ITEM_NONE;
}


/* ************************************************ */

double evalh(string model, const colvec &theta,
	     const mat &y,const mat &eta, int clustersize,
	     const SEXP &modelpar) {  
  Language call(model,
		Named("data",wrap(y)),
		Named("eta",wrap(eta)),
		Named("clustersize",wrap(clustersize)),
		Named("modelpar",modelpar));
  double h = Rcpp::as<double>(call.eval());
  return(h);
}

/* }}} */

/* ************************************************ */
/* MODELS:                                          */
/* ************************************************ */

/* {{{ weibullmm */

double weibullmm(string model, const colvec &theta,
		  const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  unsigned nx = static_cast<unsigned>(REAL(getListElement(modelpar, "nx"))[0]);
  unsigned nz = static_cast<unsigned>(REAL(getListElement(modelpar, "nlatent"))[0]);
  unsigned npar = theta.n_rows;

  double T = y(0); 
  unsigned Delta = static_cast<unsigned>(y(1)); 
  colvec X = trans(y.cols(2,nx+1));
  colvec Z = trans(y.cols(nx+2,nx+2+nz-1));
  double lambda = exp(theta(0));
  double p = exp(theta(1));
  colvec beta = theta.rows(2,nx+1);
  //  mat L = chol(reshape(theta.rows(nx+2,npar-1),nz,nz));    
  mat W = inv(diagmat(theta.rows(nx+2,npar-1)));

  double nu;
  if (nx>0)    
    nu = as_scalar(trans(X)*beta + eta*Z);
  else 
    nu = as_scalar(eta*Z);
  double part1 = Delta*log(lambda*p) + Delta*(p-1)*log(lambda*T) + 
    Delta*nu - pow(lambda*T,p)*exp(nu);
  double part2 = as_scalar(-0.5*eta*W*trans(eta))/clustersize ;
  /*  
  cerr << "T=" << T << endl;
  cerr << "Delta=" << Delta << endl;
  cerr << "X=" << endl << X;
  cerr << "Z=" << endl << Z;
  cerr << "lambda=" << lambda << endl;
  cerr << "p=" << p << endl;
  cerr << "beta=" << endl << beta;
  cerr << "W=" << endl << W;
  cerr << "part1="<<part1<<endl;
  cerr << "part2="<<part2<<endl;
  */  
  return(part1+part2);
}


double weibullmm2(string model, const colvec &theta,
		 const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  unsigned nx = static_cast<unsigned>(REAL(getListElement(modelpar, "nx"))[0]);
  unsigned nz = static_cast<unsigned>(REAL(getListElement(modelpar, "nlatent"))[0]);
  unsigned npar = theta.n_rows;
  double T = y(0); unsigned Delta = static_cast<unsigned>(y(1)); 
  colvec X = trans(y.cols(2,nx+1));
  colvec Z = trans(y.cols(nx+2,nx+2+nz-1));
  
  double lambda = exp(theta(0));
  double p = exp(theta(1));
  colvec beta = theta.rows(2,nx+1);
  //  mat L = chol(reshape(theta.rows(nx+2,npar-1),nz,nz));  
  mat L = diagmat(theta.rows(nx+2,npar-1));
  mat ceta = L*trans(eta);
  double nu = as_scalar(trans(X)*beta + trans(Z)*ceta);
  double part1 = Delta*log(lambda*p) + Delta*(p-1)*log(lambda*T) + 
    Delta*nu - pow(lambda*T,p)*exp(nu);
  double part2 = as_scalar(-0.5*eta*trans(eta))/clustersize ;
  return(part1+part2);
}


double weibullmm3(string model, const colvec &theta,
		  const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  unsigned nx = static_cast<unsigned>(REAL(getListElement(modelpar, "nx"))[0]);
  unsigned nz = static_cast<unsigned>(REAL(getListElement(modelpar, "nlatent"))[0]);
  unsigned npar = theta.n_rows;

  double T = y(0); 
  unsigned Delta = static_cast<unsigned>(y(1)); 
  colvec X = trans(y.cols(2,nx+1));
  colvec Z = trans(y.cols(nx+2,nx+2+nz-1));
  double lambda = exp(theta(0));
  double p = exp(theta(1));
  colvec beta = theta.rows(2,nx+1);
  //  mat L = chol(reshape(theta.rows(nx+2,npar-1),nz,nz));    
  mat W = inv(diagmat(theta.rows(nx+2,npar-1)));

  double nu = as_scalar(trans(X)*beta + eta*Z);
  double part1 = Delta*log(lambda*p) + Delta*(p-1)*log(lambda*T) + 
    Delta*nu - pow(lambda*T,p)*exp(nu);
  double part2 = as_scalar(-0.5*eta*W*trans(eta))/clustersize ;
  /*
  cerr << "T=" << T << endl;
  cerr << "Delta=" << Delta << endl;
  cerr << "X=" << endl << X;
  cerr << "Z=" << endl << Z;
  cerr << "lambda=" << lambda << endl;
  cerr << "p=" << p << endl;
  cerr << "beta=" << endl << beta;
  cerr << "W=" << endl << W;
  cerr << "part1="<<part1<<endl;
  cerr << "part2="<<part2<<endl;
  */
  return(part1+part2);
}


/* }}} */

/* {{{ weibull */

double weibull(string model, const colvec &theta,
	      const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  //  cerr << "theta=" << theta << endl;
  //  cerr << "y=" << y << endl;
  double lambda = 1/exp(theta(0));
  double p = exp(theta(1));
  double beta = theta(2);
  double sigma = theta(3);
  double T = y(0); unsigned Delta = y(1); double X = y(2);
  double myeta = eta(0);
  double nu = beta*X + sigma*myeta;
  double part1 = Delta*log(lambda*p) + Delta*(p-1)*log(lambda*T) + 
    Delta*nu - pow(lambda*T,p)*exp(nu);
  double part2 = (-0.5*myeta*myeta)/clustersize ;
  return(part1+part2);
}

double weibull2(string model, const colvec &theta,
	      const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  //  cerr << "theta=" << theta << endl;
  //  cerr << "y=" << y << endl;
  double lambda = 1/exp(theta(0));
  double p = exp(theta(1));
  double beta = theta(2);
  double sigma = theta(3);
  double T = y(0); unsigned Delta = y(1); double X = y(2);
  double myeta = eta(0);
  double nu = beta*X + myeta;
  double part1 = Delta*log(lambda*p) + Delta*(p-1)*log(lambda*T) + 
    Delta*nu - pow(lambda*T,p)*exp(nu);
  double part2 = (-0.5*myeta*myeta)/(sigma*sigma)/clustersize ;
  return(part1+part2);
}

/* }}} */

/* ************************************************ */

/* {{{ mimic */

double mimic(string model, const colvec &theta,
	     const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  unsigned ny = static_cast<unsigned>(REAL(getListElement(modelpar, "ny"))[0]);
  unsigned nx = static_cast<unsigned>(REAL(getListElement(modelpar, "npred"))[0]);
  unsigned npar = theta.n_rows;
  
  colvec H(ny+1);
  colvec mu = theta.rows(0,ny-1);
  colvec lambda = theta.rows(ny,2*ny-2);
  colvec beta = theta.rows(2*ny-1,2*ny-2+nx);
  mat W = inv(diagmat(theta.rows(2*ny+nx-1,npar-1)));  

  /*
  cerr << "mu=" << endl << mu;
  cerr << "lambda=" << endl << lambda;
  cerr << "beta=" << endl << beta;
  cerr << "W=" << endl << W;
  */
  /*  H(0) = y(0) - theta(0) - eta(0);
  for (unsigned i=1; i<ny; i++) {
    H(i) = y(i) - theta(i) - theta((ny-1)+i)*eta(0);
  }
  double zeta0 = eta(0);
  for (unsigned i=0; i<nx; i++) {
    zeta0 -= theta(npar-nx+i)*y(ny+i);
  }
  */
  H(0) = y(0) - mu(0) - eta(0);
  for (unsigned i=1; i<ny; i++) {
    H(i) = y(i) - mu(i) - lambda(i-1)*eta(0);
  }
  double zeta0 = eta(0);
  for (unsigned i=0; i<nx; i++) {
    zeta0 -= beta(i)*y(ny+i);
  }
  H(ny) = zeta0;    
  /* cerr << "H=" << endl << H; */
  /* cerr << "W=" << endl << W; */
  /* cerr << "npar="<<npar<<endl; */
  /* cerr << "nx="<<nx<<endl; */
  /* cerr << "ny="<<ny<<endl; */ 
  double val = as_scalar((trans(H)*W*(H)));
  return(-0.5*val);
}

double mimic2(string model, const colvec &theta,
	     const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  unsigned ny = static_cast<unsigned>(REAL(getListElement(modelpar, "ny"))[0]);
  unsigned nx = static_cast<unsigned>(REAL(getListElement(modelpar, "npred"))[0]);
  //  unsigned npar = theta.n_rows;
  
  colvec H(ny+1);
  colvec mu = theta.rows(0,ny-1);
  colvec lambda = theta.rows(ny,2*ny-1);
  colvec beta = theta.rows(2*ny,2*ny-1+nx);
  mat W(ny+1,ny+1); W.fill(0);
  //inv(diagmat(theta.rows(2*ny+nx-1,npar-1)));  
  W(ny,ny) = 1;
  for (unsigned i=0; i<ny; i++) {
    W(i,i) = 1.0/theta(2*ny+nx+i);
  }

  cerr << "mu=" << endl << mu;
  cerr << "lambda=" << endl << lambda;
  cerr << "beta=" << endl << beta;
  cerr << "W=" << endl << W;
    
  for (unsigned i=0; i<ny; i++) {
    H(i) = y(i) - mu(i) - lambda(i)*eta(0);
  }
  double zeta0 = eta(0);
  for (unsigned i=0; i<nx; i++) {
    zeta0 -= beta(i)*y(ny+i);
  }
  H(ny) = zeta0;    
  /* cerr << "H=" << endl << H; */
  /* cerr << "W=" << endl << W; */
  /* cerr << "npar="<<npar<<endl; */
  /* cerr << "nx="<<nx<<endl; */
  /* cerr << "ny="<<ny<<endl; */ 
  double val = as_scalar((trans(H)*W*(H)));
  return(-0.5*val);
}

/* }}} */

/* ************************************************ */

/* {{{ nsem1 */

double nsem1(string model, const colvec &theta,
	     const mat &y, const mat &eta, int clustersize, const SEXP &modelpar) {
  
  unsigned ny1 = static_cast<unsigned>(REAL(getListElement(modelpar, "ny1"))[0]);
  unsigned ny2 = static_cast<unsigned>(REAL(getListElement(modelpar, "ny2"))[0]);
  unsigned nx = static_cast<unsigned>(REAL(getListElement(modelpar, "npred"))[0]);
  unsigned npar = theta.n_rows;
  
  colvec H(ny1+ny2+3);
  colvec mu1 = theta.rows(0,ny1-1); 
  colvec mu2 = theta.rows(ny1,ny1+ny2-1); 
  colvec lambda1 = theta.rows(ny1+ny2,2*ny1+ny2-2);  
  colvec lambda2 = theta.rows(2*ny1+ny2-1,2*(ny1+ny2)-3);
  colvec beta0 = theta.rows(2*(ny1+ny2-1), 2*(ny1+ny2)+nx-3);
  double gamma1 = theta(2*(ny1+ny2-1)+nx);
  double gamma2 = theta(2*(ny1+ny2)-1+nx);
  mat S = diagmat(theta.rows(2*(ny1+ny2)+nx,npar-1));
  mat W = inv(S);  
  /*
  cerr << "mu1=\n" << mu1;
  cerr << "mu2=\n" << mu2;
  cerr << "lambda1=\n" << lambda1;
  cerr << "beta0=" << endl << beta0;
  cerr << "gamma1=" << gamma1 << endl;
  cerr << "gamma2=" << gamma2 <<endl;
  cerr << S;
  */

  H(0) = y(0) - mu1(0) - eta(1);
  //  H(0) = y(0) - eta(1);
  for (unsigned i=1; i<ny1; i++) {
    H(i) = y(i) - mu1(i) - lambda1(i-1)*eta(1);
  }
  H(ny1) = y(ny1) - mu2(0) - eta(2);
  //  H(ny1) = y(ny1) - eta(2);
  for (unsigned i=1; i<ny2; i++) {
    H(ny1+i) = y(ny1+i) - mu2(i) - lambda2(i-1)*eta(2);
  }

  double zeta0 = eta(0);
  for (unsigned i=0; i<nx; i++) {
    zeta0 -= beta0(i)*y(ny1+ny2+i);
  }
  H(ny1+ny2) = zeta0;    
  H(ny1+ny2+1) = eta(1)-eta(0);  
  //  H(ny1+ny2+1) = eta(1)-eta(0) - mu1(0);  
  H(ny1+ny2+2) = eta(2)-gamma1*(eta(0)) - gamma2*(eta(0)*eta(0));   
  //  H(ny1+ny2+2) = eta(2)-gamma1*(eta(0)) - gamma2*(eta(0)*eta(0)) - mu2(0);   
  double val = -0.5*as_scalar((trans(H)*W*(H))); 
  double val0; double sign0; 
  //  cerr << "S=" << S << endl;
  log_det(val0,sign0,S);
  //  cerr << "val0=" << val0 << endl;
  return(-0.5*val0+val);
}

/* }}} */

/* ************************************************ */
/* ************************************************ */
/* ************************************************ */

/* {{{ Eval */
RcppExport SEXP Eval(const SEXP modelpar, 
		     const SEXP eta,
		     const SEXP data,
		     const SEXP cluster);

SEXP Eval(const SEXP modelpar, 
	   const SEXP eta,
	   const SEXP data,
	   const SEXP cluster) {
  
  Rcpp::NumericVector theta(getListElement(modelpar, "theta"));
  colvec Theta(theta.begin(), theta.size(), 1, false); // Avoid copying

  //  List ModelPar(modelpar);
  string model = CHAR(STRING_ELT(getListElement(modelpar,"model"), 0));

  DesignFunPt modelPt = DensFixed(model);
  Rcpp::NumericMatrix D(data);
  int nobs = D.nrow(), k = D.ncol();
  mat Data(D.begin(), nobs, k, false); // Avoid copying
  Rcpp::NumericMatrix Eta(eta);  
  mat etas = mat(Eta.begin(), Eta.nrow(), Eta.ncol(), false);
  Rcpp::NumericVector C(cluster);
  colvec fCluster(C.begin(),C.size(),false);
  ucolvec Cluster = conv_to<ucolvec>::from(fCluster-1);
  int ncluster = (*max_element(Cluster.begin(),Cluster.end()))+1;
  Col<unsigned> ClusterSize(ncluster); ClusterSize.fill(0);
  for (int i=0; i<nobs; i++) {
    ClusterSize(Cluster(i))++;
  }  
  vector<double> val(nobs);
  for (int i=0; i<nobs; i++) {
    val[i] = modelPt(model, Theta,
		     Data.row(i), etas.row(Cluster(i)), ClusterSize(Cluster(i)), modelpar);    
  }
  List res;
  res["x"] = val;
  return(res);
}

/* }}} */

/* ************************************************ */

#endif /* MODELS_H */


