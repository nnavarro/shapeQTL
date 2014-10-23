//
//  shapeQTL_lm.c
//  This is a local copy of Cdqrls of the stat package
//  It will ensure portability of shapeQTL code as with R version 3.1.1
//  a 4th argument 'chk' appears in this C wrapper
//
//  Created by Nicolas Navarro on 23/10/14.
//
//

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif

/* Declarations for .Call entry points */
SEXP CdqrlsShapeQTL(SEXP x, SEXP y, SEXP tol, SEXP chk);


/* This is a local copy of Cdqrls use in lm.fit but not export. This will ensure portability of shapeQTL code
    With R version 3.1.1 a 4th argument chk appears in this C function
*/

// ---------------------------
 /* A wrapper to replace
 
 z <- .Fortran("dqrls",
 qr = x, n = n, p = p,
 y = y, ny = ny,
 tol = as.double(tol),
 coefficients = mat.or.vec(p, ny),
 residuals = y, effects = y, rank = integer(1L),
 pivot = 1L:p, qraux = double(p), work = double(2*p),
 PACKAGE="base")
 
 with less allocation.
 */

/*
 #define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
static const R_CMethodDef CEntries[]  = {
    {NULL, NULL, 0}
};
static const R_CallMethodDef CallEntries[] = {
    CALLDEF(CdqrlsShapeQTL, 4),
    {NULL, NULL, 0}
};
static const R_FortranMethodDef FortEntries[] = {
    {NULL, NULL, 0}
};
static const R_ExternalMethodDef ExtEntries[] = {
    {NULL, NULL, 0}
};

void attribute_visible R_init_shapeQTL(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortEntries, ExtEntries);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);

}
*/
SEXP CdqrlsShapeQTL(SEXP x, SEXP y, SEXP tol, SEXP chk)
{
    SEXP ans;
    SEXP qr, coefficients, residuals, effects, pivot, qraux;
    int n, ny = 0, p, rank, nprotect = 4, pivoted = 0;
    double rtol = asReal(tol), *work;
    Rboolean check = asLogical(chk);
    
    ans = getAttrib(x, R_DimSymbol);
    if(check && length(ans) != 2) error(_("'x' is not a matrix"));
    int *dims = INTEGER(ans);
    n = dims[0]; p = dims[1];
    if(n) ny = (int)(XLENGTH(y)/n); /* y :  n x ny, or an n - vector */
    if(check && n * ny != XLENGTH(y))
        error(_("dimensions of 'x' (%d,%d) and 'y' (%d) do not match"),
              n,p, XLENGTH(y));
    
    /* These lose attributes, so do after we have extracted dims */
    if (TYPEOF(x) != REALSXP) {
        PROTECT(x = coerceVector(x, REALSXP));
        nprotect++;
    }
    if (TYPEOF(y) != REALSXP) {
        PROTECT(y = coerceVector(y, REALSXP));
        nprotect++;
    }
    
    double *rptr = REAL(x);
    for (R_xlen_t i = 0 ; i < XLENGTH(x) ; i++)
        if(!R_FINITE(rptr[i])) error(_("NA/NaN/Inf in '%s'"), "x");
    
    rptr = REAL(y);
    for (R_xlen_t i = 0 ; i < XLENGTH(y) ; i++)
        if(!R_FINITE(rptr[i])) error(_("NA/NaN/Inf in '%s'"), "y");
    
    const char *ansNms[] = {"qr", "coefficients", "residuals", "effects",
        "rank", "pivot", "qraux", "tol", "pivoted", ""};
    PROTECT(ans = mkNamed(VECSXP, ansNms));
    SET_VECTOR_ELT(ans, 0, qr = duplicate(x));
    coefficients = (ny > 1) ? allocMatrix(REALSXP, p, ny) : allocVector(REALSXP, p);
    PROTECT(coefficients);
    SET_VECTOR_ELT(ans, 1, coefficients);
    SET_VECTOR_ELT(ans, 2, residuals = duplicate(y));
    SET_VECTOR_ELT(ans, 3, effects = duplicate(y));
    PROTECT(pivot = allocVector(INTSXP, p));
    int *ip = INTEGER(pivot);
    for(int i = 0; i < p; i++) ip[i] = i+1;
    SET_VECTOR_ELT(ans, 5, pivot);
    PROTECT(qraux = allocVector(REALSXP, p));
    SET_VECTOR_ELT(ans, 6, qraux);
    SET_VECTOR_ELT(ans, 7, tol);
    
    work = (double *) R_alloc(2 * p, sizeof(double));
    F77_CALL(dqrls)(REAL(qr), &n, &p, REAL(y), &ny, &rtol,
                    REAL(coefficients), REAL(residuals), REAL(effects),
                    &rank, INTEGER(pivot), REAL(qraux), work);
    SET_VECTOR_ELT(ans, 4, ScalarInteger(rank));
    for(int i = 0; i < p; i++)
        if(ip[i] != i+1) { pivoted = 1; break; }
    SET_VECTOR_ELT(ans, 8, ScalarLogical(pivoted));
    UNPROTECT(nprotect);
    
    return ans;
}


