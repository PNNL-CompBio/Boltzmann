#ifndef _DGETRS_H_
#define _DGETRS_H_ 1
extern void dgetrs(char *trans, int *n, int *nrhs, double *a, int *lda, 
		   int *ipivot, double *b, int *ldb, int *info);
#endif
