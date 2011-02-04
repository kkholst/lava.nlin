#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <iostream> 
#include <cmath>
#include <algorithm>
#include <functional>
#include <cstring>

using namespace Rcpp;
using namespace std;
using namespace arma;
#include <vector>

struct hObj
{
  double hSh;
  colvec h;
  colvec grad;
  mat hess;
};

SEXP getListElement(SEXP list, const char *str);

template <class T>
string numStr(T x) {
  ostringstream nmbstr; nmbstr << x;
  string ss = nmbstr.str();    
  return ss;
}

inline double logit(double p) { return(log(p/(1.0-p))); }
//inline double tigol(double z) { return(1.0/(1+1.0/exp(z))); } 
inline double tigol(double z) { return(1.0/(1.0+exp(-z))); } 
inline double logtigol(double z) { return(-log(1.0+exp(-z))); } 
inline double logoneminustigol(double z) { return(-z + logtigol(z)); } 


/*

template <typename T>
inline Matrix<T> reverse(const Matrix<T> &m) {
  Matrix<T> res = m;
  reverse(res.begin(), res.end());
  return(res);
}

template <typename T>
inline Matrix<T> sapply(const Matrix<T> &m, T (*fun)(T)){
  Matrix<T> res = m;
  for (unsigned i=0; i<m.rows(); i++) {
    for (unsigned j=0; j<m.cols(); j++) {
      res(i,j) = (*fun)(m(i,j));      
    }
  } 
  return(res);  
}


template <typename T>
inline Matrix<T> apply2(const Matrix<T> &m, unsigned doRow, Matrix<T> (*fun)(const Matrix<T>&)){
  if (doRow==1) {
    Matrix<T> r0 = (*fun)(m(0,_));
    return(r0);  
  } else {
    Matrix<T> r0 = (*fun)(m(_,0));
    return(r0);  
  }
}

//    cerr << apply(X,1, scythe::sum) << endl;
//    cerr << apply(X,2, scythe::sum) << endl;
template <typename T>
inline  Matrix<T> apply(const Matrix<T> &m, unsigned doRow,T (*fun)(const Matrix<T>&)){
  if(doRow==1){
    Matrix<T> res(m.rows(),1);
    for (unsigned i=0; i<m.rows(); i++) {
      Matrix<T> r = m(i,_);
      res[i] = (*fun)(r);      
      //      res(i,_) = fun(r);
    } 
    return(res);  
  } else {
    Matrix<T> res(m.cols(),1);
    for (unsigned i=0; i<m.cols(); i++) {
      Matrix<T> r = m(_,i);
      //      res(i,_) = fun(r);
      res[i] = (*fun)(r);
    }
    return(res);   
  }
}

template <typename T, typename FUNCTOR>
inline Matrix<T> apply(const Matrix<T> &m, unsigned doRow, FUNCTOR fun){
  //Matrix<T> apply(const Matrix<T> &m, unsigned doRow, T (*fun)(const Matrix<T>&)){
  if(doRow==1){
    Matrix<T> res(m.rows(),1);
    for (unsigned i=0; i<m.rows(); i++) {
      Matrix<T> r = m(i,_);
      //      res[i] = (*fun)(r);      
      res[i] = fun(r);
    } 
    return(res);  
  } else {
    Matrix<T> res(m.cols(),1);
    for (unsigned i=0; i<m.cols(); i++) {
      Matrix<T> r = m(_,i);
      res[i] = fun(r);
    }
    return(res);   
  }
}

 
template <typename T>
inline scythe::Matrix<T> cumsum(const scythe::Matrix<T> &M) {
  unsigned n = M.rows();  
  if (M.rows()==1)
    n = M.cols();    
  scythe::Matrix<T> res(M.rows(), M.cols());  
  res(0,scythe::_) = M(0,scythe::_);
  for (unsigned i=1; i<n; i++) {
    res(i,scythe::_) = res(i-1,scythe::_)+M(i,scythe::_);
  }
  return(res);
}


template <typename T>
inline std::string Rout2(scythe::Matrix<T> const& M) {
  std::ostringstream out;
  
  for(unsigned r = 0; r < M.rows(); ++r) {
    out << std::endl << "[" << r+1 << "] ";
    //    out << "";      
    for(unsigned int s = 0; s < M.cols(); ++s) {
      if(s > 0) out << ", ";
      //      out << (*this)(r,s);
      out << M(r,s);
    }
  }
  out << std::endl;
  return std::string(out.str());
}

template <typename T>
inline std::string Rout(scythe::Matrix<T> const& M) {
  std::ostringstream out;
  out << "rbind("; //<< std::endl;
  for(unsigned r = 0; r < M.rows(); ++r) {
    if(r > 0) out << "), " << std::endl;
    out << "c(";      
    for(unsigned int s = 0; s < M.cols(); ++s) {
      if(s > 0) out << ", ";
      //      out << (*this)(r,s);
      out << M(r,s);
    }
  }
  out << "))" << std::endl;
  return std::string(out.str());
}


inline string Rout(const SEXP X) {
  std::ostringstream out;
  SEXP dX;
  PROTECT(dX = getAttrib(X, R_DimSymbol));
  if (dX!=R_NilValue) {
    unsigned nrow= INTEGER(dX)[0], ncol=INTEGER(dX)[1];  
    out << "cbind(";
    for (unsigned i=0; i<ncol; i++) {
      out << "c(";
      for (unsigned j=0; j<nrow; j++) {
	unsigned pos = j + i*nrow;
	out << REAL(X)[pos];
	if (j<(nrow-1))
	  out << ", ";
      }
      out << ")";
    }
    out<< ")" << endl;
  } else {
    out << "c(";
    for (unsigned j=0; j<(length(X)); j++) {
      out << REAL(X)[j];
      if (j<(length(X)-1))
	out << ", ";
    }
    out << ")" << endl;
  }  
  UNPROTECT(1);
  return std::string(out.str());
}

*/

template <typename T>
inline std::string Rout(std::vector<T> const& x) {
  std::ostringstream out;
  out << "c(";
  for(unsigned r = 0; r < x.size(); ++r) {
    if(r > 0) out << ", ";
    //      out << (*this)(r,s);
    out << x[r];
  }
  out << ")" << std::endl;
  return std::string(out.str());
}

template <typename T>
inline std::string Rout(Mat<T> const& M) {
  std::ostringstream out;
  out << "rbind("; //<< std::endl;
  for(unsigned r = 0; r < M.n_rows; ++r) {
    if(r > 0) out << "), " << std::endl;
    out << "c(";      
    for(unsigned int s = 0; s < M.n_cols; ++s) {
      if(s > 0) out << ", ";
      out << M(r,s);
    }
  }
  out << "))" << std::endl;
  return std::string(out.str());
}



#endif /* UTILS_H */
