//# Taken from the now archived R package POT
//#Descripton of the last available version:
//#Package: POT
//#Version: 1.1-3
//#Date: 2012-10-30
//#Title: Generalized Pareto Distribution and Peaks Over Threshold
//#Author: Mathieu Ribatet <mathieu.ribatet@math.univ-montp2.fr>
//#    Maintainer: Mathieu Ribatet <mathieu.ribatet@math.univ-montp2.fr>
//#    Depends: R (>= 1.8.0)
//#Description: Some functions useful to perform a Peak Over Threshold
//#analysis in univariate and bivariate cases. A user's guide is
//#    available.
//#    License: GPL (>= 2)
//#    URL: http://r-forge.r-project.org/projects/pot/
//#    Repository: CRAN
//#    Repository/R-Forge/Project: pot
//#    Repository/R-Forge/Revision: 492
//#    Repository/R-Forge/DateTimeStamp: 2012-10-30 14:21:03
//#    Date/Publication: 2012-11-06 09:49:26
//#    Packaged: 2012-10-30 15:22:42 UTC; rforge
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()

//From POT.c
void gpdlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns);
void pplik(double *data, int *n, double *loc, double *scale,
	   double *shape, double *thresh, double *noy, double *dns);
void samlmu(double *x, int *nmom, int *n, double *lmom);

