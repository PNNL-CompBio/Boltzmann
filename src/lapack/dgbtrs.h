#ifndef _DGBTRS_H_
#define _DGBTRS_H_ 1
extern void dgbtrs_(char *transp, int *np, int *klp, int *kup, int *nrhsp,
		    double *abp, int *ldabp, int *ipivp, double *bp, int *ldbp, 
		    int *infop);
#endif
