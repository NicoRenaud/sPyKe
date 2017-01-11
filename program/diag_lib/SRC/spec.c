#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#include <cblas.h>
#include "spec.h"


///////////////////////////////////////////////
// Spec pencil total routine
///////////////////////////////////////////////
void spec_pencil_zinger(double *VAL_PRP, double *VECT_PRP,double *H, double *S, int nb_orb)
{
	// temporary eigenvalues
 	double *val_prp;
	val_prp = calloc(nb_orb,sizeof(double));


	//diagonalize with cblas_dggev
	//printf("\n   === start diagonalization of H using cblas_dggev\n");
	spec_pencil(VECT_PRP, val_prp, H, S, nb_orb);
	
	// copy the eigenvalues for sorting
	cblas_dcopy(nb_orb,val_prp,1,VAL_PRP,1);
	
	// sort the eigenvalues/eigenvectors
	qsort(VAL_PRP, nb_orb, sizeof(double), compare_doubles);
	reorder_vect_prp(VECT_PRP,VAL_PRP,val_prp,nb_orb);
	
	// free the memory
	free(val_prp);
	
}


/////////////////////////////////////////////////////////
//               Diagonalisation                       //
/////////////////////////////////////////////////////////
void spec_pencil(double *VECT_PRP, double *VAL_PRP, double *H, double *S, int n) 
{
  int i;
  int nn = n, one = 1, lwork, info;
  double tmp, *work; 
  double *VR, *VL, alphar[n], alphai[n], beta[n]; 

  // memory alloc
  VR = calloc(n*n,sizeof(double));
  VL = calloc(n*n,sizeof(double));
  
  
  // precomppute the optimal lwork
  lwork = -1;
  dggev_("N", "N", &nn, H, &nn, S, &nn, 
	  alphar, alphai, beta, VL, &one, VR, &one, 
	  &tmp, &lwork, &info); 
  
  //real computation
  lwork = (int) tmp; 
  work = (double *) malloc(sizeof(double)*lwork); 

  
  dggev_("V", "V", &nn, H, &nn, S, &nn, 
	alphar, alphai, beta, VL, &nn, VR, &nn, 
	work, &lwork, &info); 

  if(info == 3)
    printf("warning : the lapack routine dggev_ failed \n"); 
  
  // store the eigenvalues
  for(i=0;i<n;i++)
   VAL_PRP[i] = (alphar[i])/beta[i];
  
  //store the eigenvalues
  cblas_dcopy(n*n,VR,1,VECT_PRP,1);
  
  // free memory
   free(work);
   free(VL);
   free(VR);
 
}

////////////////////////////////////////////////////////
//            Reorder Eigenvector                     //
////////////////////////////////////////////////////////
void reorder_vect_prp(double *VECT_PRP,double *VAL_PRP,double *val_prp,int nb_orb)
{
   int i,j,k;
   double *vect_prp_cpy;
   vect_prp_cpy = calloc(nb_orb*nb_orb,sizeof(double));
   cblas_dcopy(nb_orb*nb_orb,VECT_PRP,1,vect_prp_cpy,1);
   
   
   // for all the eigenavalues
   for(i=0;i<nb_orb;i++)
   {
     // we find the good one
     for(j=0;j<nb_orb;j++)
     {
	if(VAL_PRP[i]==val_prp[j])
	{
	  for(k=0;k<nb_orb;k++)
	    VECT_PRP[i*nb_orb+k] = vect_prp_cpy[j*nb_orb+k];
	}
       
     }
  }
  
  free(vect_prp_cpy);
}

////////////////////////////////////////////////////////
//             Sorting function qsort                 //
////////////////////////////////////////////////////////
int compare_doubles (const void *X, const void *Y)
{
       double x = *((double *)X);
       double y = *((double *)Y);

       if (x > y)
       {
               return 1;
       }
       else
       {
               if (x < y)
               {
                       return -1;
               }
               else
               {
                       return 0;
               }
       }
}