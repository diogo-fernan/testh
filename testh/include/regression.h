#ifndef __TESTH_REGRESSION__
#define __TESTH_REGRESSION__

typedef struct {
	double m;
	double b;
} reg_Linear;

double reg_CoefficientOfDetermination (double *y, double *x, int n);

reg_Linear* reg_LeastSquareMeans (double *y, double *x, int n, int print);



#endif