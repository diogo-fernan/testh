#ifndef __TESTH_DISTRIBUTIONS__
#define __TESTH_DISTRIBUTIONS__

#include "proc.h"

enum TestHDistribution { 
	TestH_Unif		= -0x1,
	TestH_Gauss		= 0x0,
	TestH_Pareto	= 0x1,
	TestH_Exp		= 0x2
};

typedef enum TestHDistribution dist;

int dist_CheckDist (dist d);

double dist_PDFProb (dist d, double x);
double dist_CDFProb (dist d, double x);

double dist_UnifPDF (double a, double b);
double dist_GaussPDF (double x, double mu, double sigma);
double dist_ParetoPDF (double x, double xm, double alpha);
double dist_ExponentialPDF (double x, double lambda);

double dist_UnifCDF (double x, double a, double b);
double dist_GaussCDF (double x, double mu, double sigma);
double dist_ParetoCDF (double x, double xm, double alpha);
double dist_ExponentialCDF (double x, double lambda);

double dist_UnifInv (double a, double b);
double dist_GaussBoxMuller ();
double dist_GaussPolar ();
double dist_ParetoInv (double xm, double alpha);
double dist_ExponentialInv (double lambda);

void dist_Buckets (proc_Points *pt, int b, double *x, double *y, 
	double *buck, int f);
double* dist_CDFIncremental (proc_Points *pt);

double dist_F (double Rsq, int pt);

#endif