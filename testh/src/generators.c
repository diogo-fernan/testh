
#include "io.h"
#include "utilities.h"
#include "generators.h"
#include "distributions.h"
#include "rng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// gen Gaussian Inv --> H should be 0.5, uniform RNG
// proc_Process* gen_GaussianInv (double prob) {}
// would require dist_Gaussian () {}


void gen_PrintHeader (
	gen generator)
{
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF && 
		gen_CheckGen (generator) != ERR) {
		io_PrintSep ();
		char g1[50], g2[10];

		switch (generator) {
			case TestH_EG:
				strcpy (g1, "ExternGen");
				strcpy (g2, "EG");
				break;
			case TestH_RF:
				strcpy (g1, "ReadFile");
				strcpy (g2, "RF");
				break;
			case TestH_RFT:
				strcpy (g1, "ReadFileTime");
				strcpy (g2, "RFT");
				break;
			case TestH_Gauss:				
				strcpy (g1, "Gaussian");
				strcpy (g2, "Gauss");
				break;
			case TestH_AR:				
				strcpy (g1, "AggRenewal");
				strcpy (g2, "AR");
				break;
			case TestH_fBmSGA:
				strcpy (g1, "fBmSequentialGenerationAlgorithm");
				strcpy (g2, "fBm-SGA");
				break;
			case TestH_3SPG:
				strcpy (g1, "SimpleSelfSimilarProcessGenerator");
				strcpy (g2, "3SPG");
				break;
			case TestH_Hosk:
				strcpy (g1, "Hosking");
				strcpy (g2, "Hosk");
				break;

			default:
				strcpy (g1, "undefined");
				strcpy (g2, "undefined");
				break;
		}

		fprintf (stdout, "%sGENERATOR%s: %s -- %s\n",
			CGRAY_BLUE, CRESET, g1, g2);
	}
}

int gen_CheckGen (
	gen generator)
{
	int i, f;
	for (i=-3, f=0; i<=4; i++) {
		if (generator == i) {
			f = 1;
			break;
		}
	}
	if (!f)
		return ERR;
	return OK;
}


proc_Process* gen_ExternGen (N, name, gen_func, typeOfSignal)
	int N;
	char *name;
	double (*gen_func) (void);
	tosig typeOfSignal;
{
	if (N <= 0 || gen_func == NULL || 
		(typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm))
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_ExternGen");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_EG);

	int p;
	double *pts = (double*) util_MemMalloc (N * sizeof (double));
	for (p=0; p<N; p++)
		pts[p] = gen_func ();
	proc_Process *proc = proc_CreateProcess (name, pts, NULL, 
		N, TestH_fGn);
	util_MemFree (pts);
	util_TimeWr (&c);
	return proc;
}

proc_Process* gen_ReadFile (path, name, typeOfSignal)
	const char *path;
	const char *name;
	tosig typeOfSignal;
{
	if (path == NULL || name == NULL || 
		strlen (path) <= 0 || strlen (name) <= 0 ||
		(typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm))
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_ReadFile");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_RF);

	FILE *f = io_FileOpen (path, "r");
	int s = 0;
	double d, *points;
	  
	points = NULL;
	for (;;) {
		if (io_FileGetNum (f, &d) == EOF)
			break;
		points = (double*) util_MemRealloc (points, ++s * sizeof (double));
		points[s-1] = d;
		// fprintf (stdout, " %d %lf\n", s, d);
	}
	  
	io_FileClose (path, f);

	if (points == NULL)
		io_PrintErr (ERR, "invalid file %s (empty) in"
			" gen_ReadFile", path);

	proc_Process *proc = proc_CreateProcess (name, points, NULL, s, 
		typeOfSignal);

	long int mem = sizeof (points) * s;
	util_MemWr (mem);
	util_MemFree (points);
	util_TimeWr (&c);	
	return proc;
}


proc_Process* gen_ReadFileTime (path, name, typeOfSignal)
	const char *path;
	const char *name;
	tosig typeOfSignal;
{
	if (path == NULL || name == NULL || 
		strlen (path) <= 0 || strlen (name) <= 0 ||
		(typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm))
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_ReadFile");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_RFT);

	FILE *f = io_FileOpen (path, "r");
	int s1 = 0, s2 = 0;
	double d, *times, *points;
	
	times = points = NULL;
	for (;;) {
		if (io_FileGetNum (f, &d) == EOF) {
			break;
		}
		times = (double*) util_MemRealloc (times, ++s1 * sizeof (double));
		times[s1-1] = d;
		
		if (io_FileGetNum (f, &d) == EOF) {
			break;
		}
		points = (double*) util_MemRealloc (points, ++s2 * sizeof (double));
		points[s2-1] = d;
	}
	
	fclose (f);
	proc_Process *proc = NULL;

	if (points != NULL && times != NULL && s1 == s2)
		proc = proc_CreateProcess (name, points, times, s1, typeOfSignal);
	else
		io_PrintErr (ERR, "invalid file %s in"
			" gen_readFileTime", path);

	long int mem = sizeof (points) * s1 +
				sizeof (times) * s2;
	util_MemWr (mem);
	util_MemFree (points);
	util_MemFree (times);
	util_TimeWr (&c);
	return proc;
}


proc_Process* gen_Gaussian (N, typeOfSignal)
	int N;
	tosig typeOfSignal;
{
	if (N <= 0 ||
		(typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm))
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_Gaussian");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_Gauss);
	
	int p;
	double *points = (double*) malloc (N * sizeof (double));
	for (p=0; p<N; p++)
		points[p] = dist_NormalBoxMuller ();

	proc_Process *proc = proc_CreateProcess ("Gaussian", points, NULL, N, 
		typeOfSignal);
	
	if (typeOfSignal == TestH_fBm)
		proc_FractionalBrownianMotion (proc);
	
	long int mem = sizeof (points) * N;
	util_MemWr (mem);
	util_MemFree (points);
	util_TimeWr (&c);
	return proc;
}

proc_Process* gen_AggRenewal (N, h, ren, xm, typeOfSignal) 
	int N;
	double h;
	int ren;
	int xm;
	tosig typeOfSignal;
{
	if (N <= 0 || ren <= 0 || xm <= 0 || 
		(typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm) ||  
		io_CheckH (h) == ERR)
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_AggRenewal");
	
	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_AR);

	int d, s;
	int *on, *off;
	double *agg; 

	on 	= (int*)    util_MemCalloc (ren, sizeof (int));
	off	= (int*)    util_MemCalloc (ren, sizeof (int));
	agg	= (double*) util_MemCalloc (N,   sizeof (double));

	for (s=0; s<ren; s++) {
		if (rng_MT19937_genrand () > 0.5) // flip a coin
			on[s] = (int) dist_ParetoInv (xm, h);
		off[s] = (int) dist_ParetoInv (xm, h);
	}

	for (d=0; d<N; d++) {
		for (s=0; s<ren; s++) {
			if (on[s] > 0) {
				agg[d] += 1.0;
				on[s]--;
			} else { // if (off[s] > 0) 
				off[s]--;
				if (off[s] == 0) {
					on[s]  = (int) dist_ParetoInv (xm, h);
					off[s] = (int) dist_ParetoInv (xm, h);
				}
			}
		}
	}
	
	proc_Process *proc = proc_CreateProcess ("AggRenewal", agg, NULL, N, 
		typeOfSignal);
	
	proc_Normalize (proc);
	/* if (typeOfSignal == TestH_fGn)
		proc_FractionalGaussianNoise (proc); */
	if (typeOfSignal == TestH_fBm)
		proc_FractionalBrownianMotion (proc);
	
	long int mem = sizeof (agg) * N +
				sizeof (on)  * ren +
				sizeof (off) * ren;
	util_MemWr (mem);
	util_MemFree (on);
	util_MemFree (off);
	util_MemFree (agg);
	util_TimeWr (&c);
	return proc;
}

static double* gen_PersistenceProbabilities (h, scale_no)
	double h;
	int scale_no;
{
	double *prob = (double*) util_MemMalloc (scale_no * sizeof (double));
	double *expv = (double*) util_MemMalloc (scale_no * sizeof (double));
	int i;
	unsigned long long int scale = 4;

	prob[0] = pow(2, 2 * h - 2);
	expv[0] = pow(2, 2 * h - 2);

	if (TestHVerbosity > TestH_MEDIUM && TestHPrintPlain == OFF) {
		fprintf (stdout, "\n\t    Persistence Probabilities\n");
		fprintf (stdout, "\tscale\t  prob\t\t  expv\n");
		fprintf (stdout, " %10d\t%.10lf\t%.10lf\n",
			2, prob[0], expv[0]);
	}

	for (i=1; i<scale_no; i++) {
		prob[i] = (2 * pow(scale, 2 * h - 2) - pow(scale / 2, 2 * h - 2))
				/ (2 * expv[i - 1] * expv[i - 1])
				+ 0.5;

		expv[i] = expv[i - 1] * prob[i];

		if (TestHVerbosity > TestH_MEDIUM && TestHPrintPlain == OFF)
			fprintf (stdout, " %10lld\t%.10lf\t%.10lf\n",
				scale, prob[i], expv[i]);

		scale = scale * 2;
	}

	util_MemFree (expv);
	return prob;
}



proc_Process* gen_fBmSequentialGenerationAlgorithm (N, h, scale_no, typeOfSignal)
	int N;
	double h;
	int scale_no;
	tosig typeOfSignal;
{
	if (N <= 0 || !io_CheckPowerOfTwo (scale_no) ||
		(typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm) ||
		io_CheckH (h) == ERR)
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_fBmSequentialGenerationAlgorithm");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_fBmSGA);
   
	int p, i;    
	double *points, *pastPoints, *prob;

	points 		= (double*) util_MemMalloc (N * sizeof (double));
	pastPoints  = (double*) util_MemMalloc (scale_no * sizeof (double));
	prob 		= gen_PersistenceProbabilities (h, scale_no);

	int scale;
	double point, previousPoint, pastPoint;
	double proba, var;

	point = previousPoint = pastPoint = 0.0;
	proba = var = 0.0;

	for (p=0; p<N; p++) {
		// fprintf (stdout, "p %d\n", p);
		scale = 2;
		if ((p == 0)
			// 1 << scale_no = 1 * 2^scale_no
			| ((p % (1 << scale_no)) == 0)) {

			point = dist_NormalPolar ();
			// fprintf (stdout, " gauss %lf\n", point);
			for (i=0; i<scale_no; i++)
				pastPoints[i] = point;
		} else {
			int index, foundPrLevel;
			index = foundPrLevel = 0;
			while ((index < scale_no)
					& (!foundPrLevel)) {
				/* fprintf (stdout, " index %d scale %d %% %d %lf\n", 
					index, scale, p % scale, scale / 2); */
				if (p % scale == scale / 2)
					foundPrLevel = 1;
				else
					index++;
				scale = scale * 2;
			}

			proba = prob[index];
			pastPoint = pastPoints[index];
			var = sqrt (proba * (4 - proba * 4));

			point = (var * dist_NormalPolar ())
					+ (pastPoint * (2 * proba - 1));
			for (i=0; i<index; i++)
				pastPoints[i] = point;
		}

		if (typeOfSignal == TestH_fBm) {
			if (p == 1)
				point = 0.0;
			point = previousPoint + point;
			previousPoint = point;
		}
		// fprintf (stdout, " point %lf\n", point);
		points[p] = point;
	}

	proc_Process *proc = proc_CreateProcess ("fBmSequentialGenerationAlgorithm", 
		points, NULL, N, typeOfSignal);
	long int mem = sizeof (points) * N +
				sizeof (prob) * scale_no * 2 +
				sizeof (pastPoints) * scale_no;
	util_MemWr (mem);
	util_MemFree (prob);
	util_MemFree (pastPoints);
	util_MemFree (points);
	util_TimeWr (&c);
	return proc;
}


proc_Process* gen_SimpleSelfSimilarProcessGenerator (
	int 	N,
	double	h,
	int 	scale_no,
	tosig 	typeOfSignal)
{
	if ((typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm) || 
		io_CheckH (h) == ERR)
		io_PrintErr (ERR, "invalid parameters in" 
			" gen_SimpleSelfSimilarProcessGenerator");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_3SPG);

	int i, s_1, p, scale;
	double point, prob, pastPoint, var;
	double *probs, *points, *mirroredLevel, *prevPoints;
	
	prob 			= pow(2, 2 * h - 2);
	points 			= (double*) util_MemMalloc (N * sizeof (double));
	probs 			= (double*) util_MemMalloc (scale_no * sizeof(double));
	mirroredLevel 	= (double*) util_MemMalloc (scale_no * sizeof(double));
	prevPoints 		= (double*) util_MemMalloc (scale_no * sizeof(double));

	double prop1 = 1;
	double prop2 = prob;
		// Persistence Probabilities
	for(i=0; i<scale_no; i++){
		prevPoints[i] = 0;
		mirroredLevel[i] = 0;
		if (i < scale_no - 1) {
			probs[i] = sqrt (prop1 - prop2);
			if (i == scale_no - 2)
				probs[scale_no - 1] = sqrt (prop2);
			prop1 *= prob;
			prop2 *= prob;
		}
	}
	for (p=0; p<N; p++) {
		scale = 2;
		s_1 = 1;
		if ((p == 0) || ((p % (1L << scale_no)) == 0)) {
			pastPoint = 0.0;
			var = 0.0;
			point = 0.0;
			for (i = 0; i < scale_no; i++) {
				if (p == 0)
					prevPoints[i] = dist_NormalBoxMuller ();
				else {
					pastPoint = mirroredLevel[i] * prevPoints[i];
					var = sqrt (2 * prob * (1 - prob));
					prevPoints[i] = (var * dist_NormalBoxMuller ()) + 
								(pastPoint * (2 * prob - 1));
				}
				point += prevPoints[i] * probs[i];
				mirroredLevel[i] = 1;
			}
		} else {
			int index = 0;
			int foundPrLevel = 0;
			
			while ((index < scale_no) && (foundPrLevel != 1)) {
				if ((p & s_1) == (scale >> 1))
					foundPrLevel = 1;
				else
					index++;
				scale = scale << 1;
				s_1 = (s_1 << 1 ) ^ 0x01;
			}

			point -= mirroredLevel[index] * probs[index] * prevPoints[index];
			mirroredLevel[index] *= -1;
			point += mirroredLevel[index] * probs[index] * prevPoints[index];

			for (i=0; i<index; i++) {
				double pointAux = mirroredLevel[i] * prevPoints[i];
				point -= mirroredLevel[i] * probs[i] * prevPoints[i];
				var = sqrt (2 * prob * (2 - 2 * prob));
				prevPoints[i] = (var * dist_NormalBoxMuller ()) + 
							(pointAux * (2 * prob - 1));
				mirroredLevel[i] = 1;
				point += mirroredLevel[i] * probs[i] * prevPoints[i];
			}
		}

		if (typeOfSignal == TestH_fBm) {
			double pt = 0.0;
			if (p == 1)
				pt = 0.0;
			pt = pastPoint + point;
			pastPoint = pt;
		}
		points[p] = point;
	}

	proc_Process *proc = proc_CreateProcess ("SimpleSelfSimilarProcessGenerator", 
		points, NULL, N, typeOfSignal);

	long int mem = sizeof (points) * N +
				sizeof (probs) * scale_no +
				sizeof (prevPoints) * scale_no +
				sizeof (mirroredLevel) * scale_no;
	util_MemWr (mem);
	util_MemFree (probs);
	util_MemFree (prevPoints);
	util_MemFree (mirroredLevel);
	util_MemFree (points);
	util_TimeWr (&c);
	return proc;
}










double covariance(long i, double h) {
  if (i == 0)
    return 1;
  else
    return (pow(i-1,2*h)-2*pow(i,2*h)+pow(i+1,2*h))/2;
}

proc_Process* gen_Hosking (N, h, typeOfSignal)
	int N;
	double h;
	tosig typeOfSignal;
{

    if ((typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm) || 
    	io_CheckH (h) == ERR)
        io_PrintErr (ERR, "invalid parameters in" 
        	" gen_Hosking");

    clock_t c;
    util_TimeIt (&c);
    gen_PrintHeader (TestH_Hosk);

    long i, j;
    int p;
    double v, point, dPoint, dPreviousPoint = 0;
    double *phi = (double *) util_MemMalloc(N * sizeof(double));
    double *psi = (double *) util_MemMalloc(N * sizeof(double));
    double *cov = (double *) util_MemMalloc(N * sizeof(double));
    double *points=(double*) util_MemMalloc(N * sizeof(double));

    v = 1;
    phi[0] = 0;
    for (i=0; i<N; i++)
        cov[i] = covariance(i, h);
    points [0] = dist_NormalBoxMuller ();
    for(p = 1 ; p <  N; p++){
        phi[p-1] = cov[p];
        for (j=0; j<p-1; j++) {
            psi[j] = phi[j];
            phi[p-1] -= psi[j]*cov[p-j-1];
        }
        phi[p-1] /= v;
        for (j=0; j<p-1; j++) {
            phi[j] = psi[j] - phi[p-1]*psi[p-j-2];
        }
        v *= (1-phi[p-1]*phi[p-1]);

        dPoint = 0;
        for (j=0; j<p; j++) {
            dPoint += phi[j]*points[p-j-1];
        }
        dPoint += sqrt(v) * dist_NormalBoxMuller ();
        if (typeOfSignal == TestH_fBm) {
            point = dPoint + dPreviousPoint;

        }else
            point = dPoint;

        points[p] = point;
	dPreviousPoint = point;
    }

    proc_Process *proc = proc_CreateProcess ("Hosking", points, NULL, N, typeOfSignal);
    util_MemFree (phi);
    util_MemFree (psi);
    util_MemFree (cov);
    util_MemFree (points);
    util_TimeWr (&c);
    return proc;
}
