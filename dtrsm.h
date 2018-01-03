#ifndef _DTRSM_H_ 
#define _DTRSM_H_ 1
extern void dtrsm_(char *side_p, char *uplo_p, char *transa_p, char *diag_p,
		  int *m_p, int *n_p, double *alpha_p, double *a,
		  int *lda_p, double *b, int *ldb_p);
#endif
