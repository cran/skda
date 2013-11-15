
#include <math.h>
#include <float.h> /* DBL_EPSILON */

#include <Rmath.h>
#include <R_ext/Applic.h>

double pz_fmin(double ax, double bx, double (*f)(double, void *),
                  void *info, double tol);
