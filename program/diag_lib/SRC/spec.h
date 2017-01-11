
#ifndef _spec_
#define _spec_


extern int dggev_(char *jobvl, char *jobvr, int *n, double *
                  a, int *lda, double *b, int *ldb, double *alphar, 
                  double *alphai, double *beta, double *vl, int *ldvl,
                  double *vr, int *ldvr, double *work, int *lwork,
                  int *info);

void spec_pencil_zinger(double *VAL_PRP, double *VECT_PRP,double *H, double *S, int nb_orb);
void spec_pencil(double *VECT_PRP, double *VAL_PRP, double *H, double *S, int n);
void reorder_vect_prp(double *VECT_PRP,double *VAL_PRP,double *val_prp,int nb_orb);
int compare_doubles (const void *X, const void *Y);

#endif