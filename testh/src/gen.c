
#include "io.h"
#include "util.h"
#include "gen.h"
#include "dist.h"
#include "stat.h"
#include "rng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


void gen_PrintHeader (
	gen g)
{
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF && 
		gen_CheckGen (g) != ERR) {
		io_PrintSep ();
		char g1[50], g2[10];

		switch (g) {
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
			case TestH_RFTA:
				strcpy (g1, "ReadFileTimeAgg");
				strcpy (g2, "RFTA");
				break;
			case TestH_Gaussian:				
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
			case TestH_Pax:
				strcpy (g1, "Paxson");
				strcpy (g2, "Pax");
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
	gen g)
{
	int i, f;
	for (i=TestH_RFTA, f=0; i<=TestH_Pax; i++) {
		if (g == i) {
			f = 1;
			break;
		}
	}
	if (!f)
		return ERR;
	return OK;
}

#define RENEWAL 10000
#define Xm 		(double) 1.0
#define S 		16

proc_Process* gen_GenProc (
	gen 	g,
	int 	n,
	double 	h) 
{
	if (gen_CheckGen (g) == ERR || n <= 0 || io_CheckH (h) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" gen_GenProc");

	switch (g) {
		case TestH_Gaussian:
			return gen_Gaussian (n, TestH_fGn);
			break;
		case TestH_AR:
			return gen_AggRenewal (n, h, RENEWAL, Xm, TestH_fGn);
			break;
		case TestH_fBmSGA:
			return gen_fBmSequentialGenerationAlgorithm (n, h, S, TestH_fGn);
			break;
		case TestH_3SPG:
			return gen_SimpleSelfSimilarProcessGenerator (n, h, S, TestH_fGn);
			break;
		case TestH_Hosk:
			return gen_Hosking (n, h, TestH_fGn);
			break;
		case TestH_Pax:
			return gen_Paxson (n, h, TestH_fGn);
			break;
		
		default:
			return NULL;
			break;
	}
}


proc_Process* gen_ExternGen (
	int 	N,
	char 	*name,
	double 	(*gen_func) (void),
	tosig 	sig)
{
	if (N <= 0 || gen_func == NULL || 
		(sig != TestH_fGn && sig != TestH_fBm))
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_ExternGen");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_EG);

	int p;
	double *pts = (double*) util_MemMalloc (N * sizeof (double));
	for (p=0; p<N; p++)
		pts[p] = gen_func ();
	proc_Process *pr = proc_CreateProcess (name, pts, NULL, 
		N, TestH_fGn);
	util_MemFree (pts);
	util_TimeWr (&c);
	return pr;
}

proc_Process* gen_ReadFile (
	const char 	*path,
	const char 	*name,
	tosig 		sig)
{
	if (path == NULL || strlen (path) <= 0 ||
		(sig != TestH_fGn && sig != TestH_fBm))
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
	}
	  
	io_FileClose (path, f);

	if (points == NULL)
		io_PrintErr (ERR, "invalid file %s (empty) in"
			" gen_ReadFile", path);

	proc_Process *pr = proc_CreateProcess (name, points, NULL, s, sig);

	long int mem = sizeof (points) * s;
	util_MemWr (mem);
	util_MemFree (points);
	util_TimeWr (&c);	
	return pr;
}

static proc_Process* gen_ReadFileTimeAux (
	const char 	*path,
	const char 	*name,
	tosig 		sig,
	int  		mult)
{
	FILE 	*f = io_FileOpen (path, "r");
	int 	s1 = 0, s2 = 0;
	double 	d, *times, *points;
	
	times = points = NULL;
	for (;;) {
		if (io_FileGetNum (f, &d) == EOF)
			break;
		times = (double*) util_MemRealloc (times, ++s1 * sizeof (double));
		times[s1-1] = d * mult;
		
		if (io_FileGetNum (f, &d) == EOF)
			break;
		points = (double*) util_MemRealloc (points, ++s2 * sizeof (double));
		points[s2-1] = d;
	}
	fclose (f);

	proc_Process *pr = NULL;
	if (points != NULL && times != NULL && s1 == s2)
		pr = proc_CreateProcess (name, points, times, s1, sig);
	else
		io_PrintErr (ERR, "invalid file %s in"
			" gen_readFileTime", path);

	long int mem = sizeof (points) * s1 +
				sizeof (times) * s2;
	util_MemWr (mem);
	util_MemFree (points);
	util_MemFree (times);
	return pr;
}


proc_Process* gen_ReadFileTime (
	const char 	*path,
	const char 	*name,
	tosig 		sig,
	int 		mult)
{
	if (path == NULL || name == NULL || 
		strlen (path) <= 0 || strlen (name) <= 0 ||
		(sig != TestH_fGn && sig != TestH_fBm) || mult <= 0)
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_ReadFileTime");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_RFT);

	proc_Process *pr = gen_ReadFileTimeAux (path, name, sig, mult);
	
	util_TimeWr (&c);
	return pr;
}


proc_Process* gen_Gaussian (
	int 	N,
	tosig 	sig)
{
	if (N <= 0 ||
		(sig != TestH_fGn && sig != TestH_fBm))
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_Gaussian");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_Gaussian);
	
	int p;
	double *points = (double*) malloc (N * sizeof (double));
	for (p=0; p<N; p++)
		points[p] = dist_GaussBoxMuller ();

	proc_Process *pr = proc_CreateProcess ("Gaussian", points, NULL, N, sig);
	
	if (sig == TestH_fBm)
		proc_FractionalBrownianMotion (pr);
	
	long int mem = sizeof (points) * N;
	util_MemWr (mem);
	util_MemFree (points);
	util_TimeWr (&c);
	return pr;
}

proc_Process* gen_AggRenewal (
	int 	N,
	double 	h,
	int 	ren,
	int 	xm,
	tosig 	sig)
{
	if (N <= 0 || ren <= 0 || xm <= 0 || 
		(sig != TestH_fGn && sig != TestH_fBm) ||  
		io_CheckH (h) == ERR)
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_AggRenewal");
	
	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_AR);

	int 	d, s;
	int 	*on, *off;
	double 	*agg; 

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
					on[s]  = (int) dist_ParetoInv (xm, 3.0-2.0*h);
					off[s] = (int) dist_ParetoInv (xm, 3.0-2.0*h);
				}
			}
		}
	}
	
	proc_Process *pr = proc_CreateProcess ("AggRenewal", agg, NULL, N, 
		sig);
	
	proc_Normalize (pr);
	if (sig == TestH_fBm)
		proc_FractionalBrownianMotion (pr);
	
	long int mem = sizeof (agg) * N +
				sizeof (on)  * ren +
				sizeof (off) * ren;
	util_MemWr (mem);
	util_MemFree (on);
	util_MemFree (off);
	util_MemFree (agg);
	util_TimeWr (&c);
	return pr;
}

static double* gen_PersistenceProbabilities (
	double 	h,
	int 	scale_no)
{
	double *prob = (double*) util_MemMalloc (scale_no * sizeof (double));
	double *expv = (double*) util_MemMalloc (scale_no * sizeof (double));
	int i;
	unsigned long long scale = 4;

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

proc_Process* gen_fBmSequentialGenerationAlgorithm (
	int 	N,
	double 	h,
	int 	scale_no,
	tosig 	sig)
{
	if (N <= 0 || !io_CheckPowerOfTwo (scale_no) ||
		(sig != TestH_fGn && sig != TestH_fBm) ||
		io_CheckH (h) == ERR)
		  io_PrintErr (ERR, "invalid parameters in"
			  " gen_fBmSequentialGenerationAlgorithm");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_fBmSGA);
   
	int 	p, i, m;    
	double 	proba, var;
	double 	point, prevPoint, pastPoint;
	double 	*points, *pastPoints, *prob;

	points 		= (double*) util_MemMalloc (N * sizeof (double));
	pastPoints  = (double*) util_MemMalloc (scale_no * sizeof (double));
	prob 		= gen_PersistenceProbabilities (h, scale_no);

	point = prevPoint = pastPoint = 0.0;
	proba = var = 0.0;

	for (p=0; p<N; p++) {
		// fprintf (stdout, "p %d\n", p);
		m = 2;
		if ((p == 0)
			// 1 << scale_no = 1 * 2^scale_no
			| ((p % (1 << scale_no)) == 0)) {

			point = dist_GaussPolar ();
			// fprintf (stdout, " gauss %lf\n", point);
			for (i=0; i<scale_no; i++)
				pastPoints[i] = point;
		} else {
			int index, foundPrLevel;
			index = foundPrLevel = 0;
			while ((index < scale_no)
					& (!foundPrLevel)) {
				/* fprintf (stdout, " index %d m %d %% %d %lf\n", 
					index, m, p % m, m / 2); */
				if (p % m == m / 2)
					foundPrLevel = 1;
				else
					index++;
				m *= 2;
			}

			proba = prob[index];
			pastPoint = pastPoints[index];
			var = sqrt (proba * (4 - proba * 4));

			point = (var * dist_GaussPolar ())
					+ (pastPoint * (2 * proba - 1));
			for (i=0; i<index; i++)
				pastPoints[i] = point;
		}

		if (sig == TestH_fBm) {
			point 	 += prevPoint;
			prevPoint = point;
		}
		points[p] = point;
	}

	proc_Process *pr = proc_CreateProcess ("fBmSequentialGenerationAlgorithm", 
		points, NULL, N, sig);
	long int mem = sizeof (points) * N +
				sizeof (prob) * scale_no * 2 +
				sizeof (pastPoints) * scale_no;
	util_MemWr (mem);
	util_MemFree (prob);
	util_MemFree (pastPoints);
	util_MemFree (points);
	util_TimeWr (&c);
	return pr;
}

proc_Process* gen_SimpleSelfSimilarProcessGenerator (
	int 	N,
	double	h,
	int 	scale_no,
	tosig 	sig)
{
	if (N <= 0 || (sig != TestH_fGn && sig != TestH_fBm) || 
		io_CheckH (h) == ERR)
		io_PrintErr (ERR, "invalid parameters in" 
			" gen_SimpleSelfSimilarProcessGenerator");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_3SPG);

	int 	i, s_1, p, m;
	double 	point, prob, prevPoint, pastPoint, var;
	double 	*probs, *points, *mirroredLevel, *prevPoints;
	
	prob 			= pow (2, 2 * h - 2);

	points 			= (double*) util_MemMalloc (N * sizeof (double));
	probs 			= (double*) util_MemMalloc (scale_no * sizeof(double));
	mirroredLevel 	= (double*) util_MemMalloc (scale_no * sizeof(double));
	prevPoints 		= (double*) util_MemMalloc (scale_no * sizeof(double));

	prevPoint = 0.0;
	double prop1 = 1;
	double prop2 = prob;

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
		m = 2;
		s_1 = 1;
		if ((p == 0) || ((p % (1L << scale_no)) == 0)) {
			pastPoint = 0.0;
			var = 0.0;
			point = 0.0;
			for (i = 0; i < scale_no; i++) {
				if (p == 0)
					prevPoints[i] = dist_GaussBoxMuller ();
				else {
					pastPoint = mirroredLevel[i] * prevPoints[i];
					var = sqrt (2 * prob * (1 - prob));
					prevPoints[i] = (var * dist_GaussBoxMuller ()) + 
								(pastPoint * (2 * prob - 1));
				}
				point += prevPoints[i] * probs[i];
				mirroredLevel[i] = 1;
			}
		} else {
			int index = 0;
			int foundPrLevel = 0;
			
			while ((index < scale_no) && (foundPrLevel != 1)) {
				if ((p & s_1) == (m >> 1))
					foundPrLevel = 1;
				else
					index++;
				m = m << 1;
				s_1 = (s_1 << 1) ^ 0x01;
			}

			point -= mirroredLevel[index] * probs[index] * prevPoints[index];
			mirroredLevel[index] *= -1;
			point += mirroredLevel[index] * probs[index] * prevPoints[index];

			for (i=0; i<index; i++) {
				double pointAux = mirroredLevel[i] * prevPoints[i];
				point -= mirroredLevel[i] * probs[i] * prevPoints[i];
				var = sqrt (2 * prob * (2 - 2 * prob));
				prevPoints[i] = (var * dist_GaussBoxMuller ()) + 
							(pointAux * (2 * prob - 1));
				mirroredLevel[i] = 1;
				point += mirroredLevel[i] * probs[i] * prevPoints[i];
			}
		}

		if (sig == TestH_fBm) {
			point 	 += prevPoint;
			prevPoint = point;
		}
		points[p] = point;
	}

	proc_Process *pr = proc_CreateProcess ("SimpleSelfSimilarProcessGenerator", 
		points, NULL, N, sig);

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
	return pr;
}

proc_Process* gen_Hosking (
	int 	N,
	double 	h,
	tosig 	sig)
{

	if (N <= 0 || (sig != TestH_fGn && sig != TestH_fBm) || 
		io_CheckH (h) == ERR)
		io_PrintErr (ERR, "invalid parameters in" 
			" gen_Hosking");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_Hosk);

	int 	p;
	long 	i, j;
	double 	v, point, prevPoint = 0;

	double *phi = (double *) util_MemMalloc (N * sizeof(double));
	double *psi = (double *) util_MemMalloc (N * sizeof(double));
	double *cov = (double *) util_MemMalloc (N * sizeof(double));
	double *points = (double*) util_MemMalloc (N * sizeof(double));

	v = 1;
	phi[0] = 0;
	for (i=0; i<N; i++)
		cov[i] = stat_Covariance (i, h);
	points[0] = 0.0; // dist_GaussBoxMuller ();
	for(p=1; p<N; p++){
		phi[p-1] = cov[p];
		for (j=0; j<p-1; j++) {
			psi[j] = phi[j];
			phi[p-1] -= psi[j]*cov[p-j-1];
		}
		phi[p-1] /= v;
		for (j=0; j<p-1; j++)
			phi[j] = psi[j] - phi[p-1]*psi[p-j-2];
		v *= (1 - phi[p-1] * phi[p-1]);

		point = 0;
		for (j=0; j<p; j++) 
			point += phi[j] * points[p-j-1];
		point += sqrt (v) * dist_GaussBoxMuller ();

		if (sig == TestH_fBm) {
			point = prevPoint + point;
			prevPoint = point;        	
		}
		points[p] = point;
	}

	proc_Process *pr = proc_CreateProcess ("Hosking", points, NULL, N, sig);

	long int mem = sizeof (phi) * N +
				sizeof (psi) * N +
				sizeof (cov) * N +
				sizeof (points) * N;
	util_MemWr (mem);
	util_MemFree (phi);
	util_MemFree (psi);
	util_MemFree (cov);
	util_MemFree (points);
	util_TimeWr (&c);
	return pr;
}


typedef	struct	{
	unsigned int dim, max_dim;
	double *ve;
} Paxson_VEC;

static Paxson_VEC* Paxson_VGet (
	unsigned int dim)
{
	Paxson_VEC *v;
	v = (Paxson_VEC*) util_MemMalloc (sizeof (Paxson_VEC));
	v->ve = (double*) util_MemMalloc (dim * sizeof (double));
	v->dim = dim;
	v->max_dim = dim;
	return v;
}

static Paxson_VEC* Paxson_VResize (
	Paxson_VEC	 *v,
	unsigned int dim)
{
	v->ve = util_MemRealloc (v->ve, dim * sizeof (double));
	v->dim = dim;
	v->max_dim = dim;
	return v;
}

static void Paxson_FastFourierTransform (
	Paxson_VEC *x_re, 
	Paxson_VEC *x_im)
{
	if (!x_re || !x_im ||
		x_re->dim != x_im->dim)
		io_PrintErr (ERR, "invalid parameters in" 
			" Paxson_FastFourierTransform");

	int 	i, ip, j, k, li, n, length;
	double	*xr, *xi;
	double	theta;
	double 	w_re, w_im, u_re, u_im, t_re, t_im;
	double 	tmp, tmpr, tmpi;

	n = 1;
	while (x_re->dim > n)
		n *= 2;
	x_re = Paxson_VResize (x_re, n);
	x_im = Paxson_VResize (x_im, n);
	xr   = x_re->ve;
	xi   = x_im->ve;

	/* Decimation in time (DIT) algorithm */
	j = 0;
	for (i=0; i<n-1; i++) {
		if (i < j) {
			tmp   = xr[i];
			xr[i] = xr[j];
			xr[j] = tmp;
			tmp   = xi[i];
			xi[i] = xi[j];
			xi[j] = tmp;
		}
		k = n / 2;
		while (k <= j) {
			j -= k;
			k /= 2;
		}
		j += k;
	}

	/* Actual FFT */
	for (li=1; li<n; li *= 2) {
		length = 2 * li;
		theta  = PI / li;
		u_re = 1.0; u_im = 0.0;
		if (li == 1) {
			w_re = -1.0;
			w_im =  0.0;
		}
		else if (li == 2) {
			w_re =  0.0;
			w_im =  1.0;
		}
		else {
			w_re = cos (theta);
			w_im = sin (theta);
		}
		for (j=0; j<li; j++) {
			for (i=j; i<n; i += length) {
				ip = i + li;
				/* step 1 */
				t_re = xr[ip]*u_re - xi[ip]*u_im;
				t_im = xr[ip]*u_im + xi[ip]*u_re;
				/* step 2 */
				xr[ip] = xr[i] - t_re;
				xi[ip] = xi[i] - t_im;
				/* step 3 */
				xr[i] += t_re;
				xi[i] += t_im;
			}
			tmpr = u_re*w_re - u_im*w_im;
			tmpi = u_im*w_re + u_re*w_im;
			u_re = tmpr;
			u_im = tmpi;
		}
	}
}

static double Paxson_fGn_BEstAdj (
	double lambda, 
	double H) 
{
	int 	k;
	double 	d, prime, sum1, sum2, result;
	double 	a[5], b[5]; /* index 0 never used ! */

	d = -2.0 * H - 1.0;
	prime = -2.0 * H;

	for(k=1; k<5; k++) {
		a[k] = 2.0 * k * PI + lambda;
		b[k] = 2.0 * k * PI - lambda;
	}
	sum1 = 0.0;
	for(k=1; k<4; k++) {
		sum1 += pow (a[k], d);
		sum1 += pow (b[k], d);
	}
	sum2 = 0.0;
	for(k=3; k<5; k++) {
		sum2 += pow (a[k], prime);
		sum2 += pow (b[k], prime);
	}
	result = sum1 + (sum2 / (8.0 * PI * H));
	return (1.0002 - 0.000134 * lambda) *
			(result - pow(2, -7.65 * H - 7.4));
}

static void Paxson_fGn_Spectrum (
	double 	*pow_spec, 
	int 	n, 
	double 	H) 
{
	int 	 i;
	double lambda, fact1, a, b, c, g;

	/* the result of lgamma will always be positive: */
	g = lgamma (2.0 * H + 1.0);
	fact1 = 2.0 * sin (PI * H) * exp(g);

	for(i=1; i<n+1; i++) {
		lambda = (PI * i) / n;
		a = fact1 * (1.0 - cos (lambda));
		b = pow (lambda, -2.0 * H - 1.0);
		c = Paxson_fGn_BEstAdj (lambda, H);
		pow_spec[i] = a * (b + c);
	}
}

proc_Process* gen_Paxson (
	int 	N,
	double 	h,
	tosig 	sig)
{

	if (N <= 0 || (sig != TestH_fGn && sig != TestH_fBm) || 
		io_CheckH (h) == ERR)
	    io_PrintErr (ERR, "invalid parameters in" 
	    	" gen_Paxson");

	clock_t c;
	util_TimeIt (&c);
	gen_PrintHeader (TestH_Pax);

	int 	p;
	long 	halfn;
	double 	aux;
	Paxson_VEC *a_re, *a_im;

	halfn = N / 2;

	double *pow_spec = (double*) util_MemMalloc ((halfn+1) * sizeof (double));
	double *points 	 = (double*) util_MemMalloc (N * sizeof (double));

	/* approximate spectral density */
	Paxson_fGn_Spectrum (pow_spec, halfn, h);

	a_re = Paxson_VGet (N);
	a_im = Paxson_VGet (N);
	a_re->ve[0] = 0;
	a_im->ve[0] = 0;
	for(p=1; p<=halfn; p++) {
		aux = sqrt (pow_spec[p]);
		a_re->ve[p] = aux * dist_GaussBoxMuller ();
		a_im->ve[p] = aux * dist_GaussBoxMuller ();
	}
	for(p=halfn+1; p<N; p++) {
		a_re->ve[p] =  a_re->ve[N-p];
		a_im->ve[p] = -a_im->ve[N-p];
	}

	/* real part of Fourier transform of a_re + i a_im gives sample path */
	Paxson_FastFourierTransform (a_re, a_im);

	points[0] = a_re->ve[0];
	for(p=1; p<N; p++) {
		points[p] = a_re->ve[p];
		if (sig == TestH_fBm)
			points[p] += points[p-1];
	}

	proc_Process *pr = proc_CreateProcess ("Paxson", points, NULL, N, sig);

	long int mem = sizeof (pow_spec) * (halfn+1) +
				sizeof (points) * N;
	util_MemWr (mem);
	util_MemFree (pow_spec);
	util_MemFree (points);
	util_TimeWr (&c);
	return pr;
}








/* Random Midpoint Displacement (RMD)

Dear Bernardo,

nice to hear of your and Prof. Inácio's interest in this work. Indeed, we have not kept our old C code alive. However, Ton Dieker has a wonderful fBm page http://www2.isye.gatech.edu/~adieker3/fbm.html that includes also our algorithm - please look there.

BTW, nobody (as far as I know) has characterized the accuracy of our approximate algorithm through some metric concerning the whole covariance structure, not only the Hurst parameter. It would be nice to know whether for example m=n=10 is "practically perfect" or not. If you find some fresh idea in that direction, please send a preprint to me also!

Best regards,
Ilkka Norros */
