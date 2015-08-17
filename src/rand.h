double sTruncNorm(double bd, double mu, double var, int lower);
double TruncNorm(double lb, double ub, double mu, double var, int invcdf);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
void rWish(double **Sample, double **S, int df, int size);
void RTruncNorm(
     double * x,
		 double * lb,  /* lower bound */ 
		 double * ub,  /* upper bound */
		 double * mu,  /* mean */
		 double * var, /* variance */
		 int * invcdf  /* use inverse cdf method? */
  ); 

void RWISH(
    double * sample,
    double * s,
    int * dfPtr,
    int * sizePtr
    );

void RMVN(
    double * sample,
    double * mean,
    double * var,
    int * sizePtr
    ); 


