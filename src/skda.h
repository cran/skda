
#include <R.h>
#include <R_ext/Lapack.h>





typedef struct
{
  const int *p;
  const int *n;
  const int *nclass;
  const int *y;
  const double *dsqx;
  const double *prop;
  double *lamtmp0;
  double *lamtmp1;
  double *lam;
  double *phat;
  int *ny;
  int *nsqp; 
  int *np;
  double tau;
} skdaCDinfo;



double skdaCD_obj(double c, void *arg);


void coordescent(const int *n, const int *p, const int *nclass, const double *x,
                 const int *y, const double *tau,
                 const int *maxct, const double *prop, double *lam, double *phat);

void condprob(const int *n, const int *p, const int *nclass,
              const int *y, const double *x, const double *lam,
              const int *m, const double *newx, const double *prop,
              double *phat);

void condSMXprob(const int *n, const int *p, const int *nclass,
              const int *y,  const double *lam,
              const double *dsqx, const int *nsqp, const int *np,
              const int *ny, const double *prop,
              double *phat);

