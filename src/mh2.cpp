#include "utils.h"
#include "mh2.h"
#include "models.h"
//#include <time.h>

SEXP MH(
	SEXP data,  
	SEXP cluster,
	SEXP etainit,
	SEXP Sigma,
	SEXP modelpar,
	SEXP control
	) {   
  
  RNGScope scope;
  //  srand ( time(NULL) ); /* initialize random seed: */

  Rcpp::NumericVector theta(getListElement(modelpar, "theta"));
  colvec Theta(theta.begin(), theta.size(), 1, false); // Avoid copying
  Rcpp::NumericMatrix D(data);
  int nobs = D.nrow(), k = D.ncol();
  mat Data(D.begin(), nobs, k, false); // Avoid copying
  Rcpp::NumericMatrix V(Sigma);  

  mat VarEta(V.begin(), V.nrow(), V.ncol(), false); 

  Rcpp::NumericVector C(cluster);
  colvec fCluster(C.begin(),C.size(),false);
  ucolvec Cluster = conv_to<ucolvec>::from(fCluster-1);

  mat L = chol(VarEta);  
  int ncluster = (*max_element(Cluster.begin(),Cluster.end())) +1;
  Col<unsigned> ClusterSize(ncluster); ClusterSize.fill(0);
  for (int i=0; i<nobs; i++) {
    ClusterSize(Cluster(i))++;
  }  

  List ModelPar(modelpar);
  List MCMCopt(control);
  IntegerVector _nlatent = ModelPar["nlatent"]; unsigned nlatent = _nlatent[0];
  IntegerVector _internal = ModelPar["internal"]; bool internal = _internal[0];
  string model = CHAR(STRING_ELT(getListElement(modelpar,"model"), 0));  
  IntegerVector _m = MCMCopt["m"]; unsigned m = _m[0];
  NumericVector _stepsize = MCMCopt["stepsize"]; double stepsize = _stepsize[0];
  IntegerVector _thin = MCMCopt["thin"]; unsigned thin = _thin[0];
  IntegerVector _nsim = MCMCopt["nsim"]; unsigned nsim = _nsim[0];
  //  RcppParams ModelPar(modelpar);
  //  RcppParams MCMCopt(control);     
  // int nlatent = ModelPar.getIntValue("nlatent");
  // string model = ModelPar.getStringValue("model");
  // bool internal = ModelPar.getBoolValue("internal");
  // double stepsize = MCMCopt.getDoubleValue("stepsize");
  // int m = MCMCopt.getIntValue("m");
  // int thin = MCMCopt.getIntValue("thin");
  // int nsim = MCMCopt.getIntValue("nsim");
  // cerr << "nlatent=" << nlatent << endl;
  // cerr << "internal=" << internal << endl;  
  // cerr << "stepsize=" << stepsize << endl;
  // cerr << "m=" << m << endl;
  // cerr << "thin=" << thin << endl;
  // cerr << "nsim=" << nsim << endl;
  // cerr << "model=" << model << endl;


  DesignFunPt modelPt = evalh; /* Default: model design from R function 
				  defined by the string: model */
  if (internal) {  
    modelPt = DensFixed(model);
  }

  //  cerr << "cl=" << Rout(ClusterSize) << endl;   
  //  cerr << "model=" << model << endl;
  //  cerr << "ncluster=" << ncluster << endl;
  //  cerr << "nlatent=" << nlatent << endl;

  //  mat Init(I.begin(), I.size(), 1, false);
  //  mat Theta = Init;
  colvec eta(nlatent);
  //mat etas = zeros(ncluster,nlatent);  
  Rcpp::NumericMatrix EtaInit(etainit);  
  mat etas = mat(EtaInit.begin(), EtaInit.nrow(), EtaInit.ncol(), false);
  //  if (etainit!=NILSXP) {
  //    Rcpp::NumericMatrix EtaInit(etainit);  
  //  }
  mat oldetas = etas;    

  rowvec logf_old(ncluster); logf_old.fill(0);
  rowvec logf_prop(ncluster);  
  for (int i=0; i<nobs; i++) {
    logf_old(Cluster(i)) += 
      modelPt(model, Theta,
	    Data.row(i), oldetas.row(Cluster(i)), ClusterSize(Cluster(i)), modelpar);    
  }
  vector<int> mwhich(m);
  for (int i=0; i<m; i++) {
    mwhich[m-1-i] = nsim-(thin+1)*i;
  }
  int mcounter = 0;

  int n = ncluster*nlatent;
  mat chain(nsim, n);
  colvec accept = zeros(ncluster); 

  mat CompleteData;
  //  cerr << "VarEta=" << VarEta;
  //  cerr << "L=" << L;
  //  cerr << "stepsize=" << stepsize << endl;
  
  
  for (int iter=0; iter<nsim; iter++) {
    for (int k=0; k<ncluster; k++) {
      logf_prop(k) = 0; 
      //      colvec Z = randn(nlatent);      
      colvec Z = rnorm(nlatent,0,1);
      eta = stepsize*(L*Z);
      etas.row(k) = trans(eta);
    }
    etas = etas + oldetas;
    for (int i=0; i<nobs; i++) {
      logf_prop(Cluster(i)) += 
	modelPt(model, Theta,
	      Data.row(i), etas.row(Cluster(i)), ClusterSize(Cluster(i)), modelpar);
    } // i (nobs)

    //    colvec U = randu(ncluster);
    colvec U = runif(ncluster,0,1);
    // cerr << "---------------------" << endl;
    // cerr << "logf_prop=" << Rout(logf_prop) << endl;
    // cerr << "logf_old=" << Rout(logf_old) << endl;
    for (int k=0; k<ncluster; k++) {
      double alpha = fmin(0, (logf_prop(k)-logf_old(k)));
      //      double alpha = fmin(0, (logf_prop(k)-logf_old(k)));
      //      cerr << "alpha=" << alpha << endl;
      //      if (log(U(k))<alpha) { // accept
//      cerr << "U=" << U(k);
//      cerr << ", alpha=" << alpha << endl;
      if (log(U(k))<alpha) { // accept
	oldetas.row(k) = etas.row(k);	// New value
	logf_old(k) = logf_prop(k);
	accept(k)++;
      }
    }       
    //    cerr << "accept = " << accept/nsim << endl;
    chain.row(iter) = reshape(oldetas,1,n);

    
    if (iter==(mwhich[mcounter]-1)) {
      mat myeta(nobs,nlatent); // Last draw
      for (int i=0; i<nobs; i++) {
	myeta.row(i) = oldetas.row(Cluster(i));
      }
      if (mcounter==0) {
	CompleteData = join_rows(Data,myeta);
      } else {
	CompleteData = join_cols(CompleteData, 
				 join_rows(Data,myeta));
      }      
      mcounter++;
    }
    

  } // iter (nsim)

  /*
  mat myeta(nobs,nlatent); // Last draw
  for (int i=0; i<nobs; i++) {
	myeta.row(i) = oldetas.row(Cluster(i));
  }
  */
  List  res;
  res["accept"] = accept;
  res["chain"] = chain;
  res["eta"] = oldetas;
  res["data"] = CompleteData;
  //  res["data"] = join_rows(Data,myeta);
  return res;
}




//##################################################

SEXP FastApprox(const SEXP a,
		const SEXP t,
		const SEXP z) {
  NumericVector A(a);
  NumericVector T(t);
  NumericVector Z(z);
  vector<unsigned> idx(Z.size());
  vector<double> newT(Z.size());

  NumericVector::iterator it;  
  double lower,upper; int pos=200;
  for (int i=0; i<Z.size(); i++) {
    it = lower_bound(A.begin(), A.end(), Z[i]);
    if (it == A.begin()) { 
      pos = 0; 
      // upper = *it; // no smaller value  than val in vector
    } 
    else if (int(it-A.end())==0) {
      pos = A.size()-1;
      //lower = *(it-1); // no bigger value than val in vector
    } else {
      lower = *(it-1);
      upper = *it;
      pos = int(it- A.begin());
      if (abs(Z[i]-lower)<abs(Z[i]-upper)) {
	pos = int(it- A.begin());
      }
    }
    idx[i] = pos;
    newT[i] = T[pos];
  }
  List ans;
  ans["t"] = newT;
  ans["pos"] = idx;
  return(ans);
}

//##################################################
