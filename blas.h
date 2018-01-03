#ifndef _BLAS_H_ 
#define _BLAS_H_
extern void dscal(int *nx, double *a_p, double *x, int *inc);
extern void dcopy(int *nx, double *x, int *incx, double *y, int *incy);
extern double dnrm2(int *nx, double *x, int *inc);
extern void dgemv(int *trans, int *m_p, int *n_p, double *alpha_p, double *a, 
		  int *lda_p, double *x, int *incx_p, double *beta_p, 
		  double *y, int *incy_p);
extern void daxpy(int *n_p, double *alpha_p, double *x, int *incx_p, 
		  double *y, int *incy_p);
#endif
