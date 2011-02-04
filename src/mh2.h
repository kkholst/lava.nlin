#ifndef _MH2_H
#define _MH2_H

//#include <R.h>           //  Rprintf()
//#include <R_ext/Utils.h> //  user interrupts
//#include <Rdefines.h>
//#include <Rinternals.h>

RcppExport SEXP MH(SEXP data,
		     SEXP cluster,
		   //		   SEXP init,
		     SEXP etainit,
		     SEXP Sigma,
		     SEXP modelpar, 
		     SEXP control);


RcppExport SEXP FastApprox(const SEXP a,
			   const SEXP t,
			   const SEXP z);


#endif /* _MH2_H */
