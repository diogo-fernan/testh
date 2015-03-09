
#include "proc.h"
#include "io.h"
#include "util.h"
#include "reg.h"
#include "stat.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static char* proc_TypeOfScale (tos typeOfScale);

proc_Process* proc_CreateProcess (
	const char 	*n,
	double 		*p,
	double 		*t,
	int 		s,
	tosig 		sig)
{
	if (p == NULL || s <= 0 || 
		(sig != TestH_fGn && sig != TestH_fBm))
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateProcess");

	proc_Process *pr = (proc_Process*) util_MemMalloc (sizeof (proc_Process));
	pr->points = proc_CreatePoints (p, t, s, 1, sig);
	pr->scales = NULL;
	pr->scales_on = OFF;
	if (n != NULL) {
		int st = strlen (n) + 1;
		pr->name = (char*) util_MemMalloc (st * sizeof (char));
		memset (pr->name, '\0', st);
		memcpy (pr->name, n, st);
	} else {
		pr->name = (char*) util_MemMalloc (11 * sizeof (char));
		memset (pr->name, '\0', 11);
		pr->name = "N/A";
	}
	return pr;
}

proc_Points* proc_CreatePoints (
	double 	*p,
	double 	*t,
	int 	s,
	int 	m,
	tosig 	sig)
{
	if (p == NULL || s < 0 || m < 0)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreatePoints");

	int i;
	proc_Points *pt = (proc_Points*) util_MemMalloc (sizeof (proc_Points));
	pt->size = s;
	pt->scale = m;
	pt->typeOfSignal = sig;
	pt->points = (double*) util_MemMalloc (s * sizeof (double));
	
	for (i=0; i<s; i++)
		pt->points[i] = p[i];
	
	if (t != NULL) {
		pt->times = (double*) util_MemMalloc (s * sizeof (double));
		for (i=0; i<s; i++)
			pt->times[i] = t[i];
	} else
		pt->times = NULL;
	return pt;
}

// ... --> for manual scaling
proc_ScalesConfig* proc_CreateScalesConfig (
	tos scl,
	int min,
	int max,
	int mov,
	...)
{
	if (scl != TestH_INC && 
		scl != TestH_RAND && 
		scl != TestH_POW)
		io_PrintErr (ERR, "invalid scale type in"
			" proc_CreateScalesConfig");
	if (min <= 0 || max <= 0 || mov <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScalesConfig");

	int i;
	long long int s;
	proc_ScalesConfig *conf = (proc_ScalesConfig*) util_MemMalloc (
		sizeof (proc_ScalesConfig));

	conf->min = min;
	conf->max = max;
	conf->mov = mov;
	conf->typeOfScale = scl;
	if (scl == TestH_INC)
		conf->no = (max - min) / mov + 1;
	// else if (scl == TestH_RAND)
	else // if (scl == TestH_POW) 
		conf->no = max - min + 1;

	io_PrintSep ();        
	// proc_PrintHeaderScales (conf);
	if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF) {
		fprintf (stdout, " - %d scales, min: %d, max: %d",
			conf->no, min, max);
		if (scl == TestH_INC)
			fprintf (stdout, ", increment: %d\n", mov);
		else if (conf->typeOfScale == TestH_RAND)
			fprintf (stdout, "\n");
		else if (conf->typeOfScale == TestH_POW)
			fprintf (stdout, ", exp base: %d\n", mov);
		fprintf (stdout, "CreateScalesConfig -- %s\n", 
			proc_TypeOfScale (conf->typeOfScale));
		fprintf (stdout, "\n [");
	}

	conf->scales = (int*) util_MemCalloc (conf->no, sizeof (int));
	for (i=0, s=min; i<conf->no; i++) {
		if (scl == TestH_POW)
			s = pow (mov, min + i);

		conf->scales[i] = s;
		if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF) {
			fprintf (stdout, " %lld", s);
			if (i != conf->no - 1)
				fprintf (stdout, ",");
		}

		if (scl == TestH_INC)
			s += mov;
		// else if (scl == TestH_RAND)
		// else if (scl == TestH_MAN)
	}
	if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF)
		fprintf (stdout, "]\n");

	if (proc_CheckScalesConfig (conf) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScalesConfig");
	return conf;
}

void proc_CreateScalesTimeAux (
	proc_Process 		*pr,
	proc_ScalesConfig 	*conf,
	int 				agg)
{
	int 			i, p, t, tfirst, tlast, tstep;
	long long int 	size;
	double 			ts, res;
	int 			*scales_sam = NULL;
	double 			*scales_tot = NULL;
	double 			*scales_pos = NULL;
	proc_Scales 	*scales = NULL;

	pr->scales = (proc_Scales*) util_MemMalloc (sizeof (proc_Scales));
	scales = pr->scales;
	scales->conf = conf;
	scales->scales = (proc_Points*) util_MemMalloc (conf->no * 
		sizeof (proc_Points));
	
	scales_sam = (int*)    util_MemCalloc (conf->no, sizeof (int));
	scales_tot = (double*) util_MemCalloc (conf->no, sizeof (double));
	scales_pos = (double*) util_MemCalloc (conf->no, sizeof (double));

	for (i=0; i<conf->no; i++) {
		scales->scales[i].size = 0;
		scales->scales[i].scale = conf->scales[i];
		scales->scales[i].typeOfSignal = TestH_fGn;
		scales->scales[i].points = NULL;
		scales->scales[i].times = NULL;

		scales_pos[i] = floor (pr->points->times[0]);
	}

	tfirst = floor (pr->points->times[0]);
	tlast  = ceil  (pr->points->times[pr->points->size-1]);
	for (i=0; i<conf->no; i++) {
		tstep = conf->scales[i];
		p = 0;
		for (t=tfirst; ; t+=tstep) {
			if (t > tlast)
				break;
			for (;;) {
				if (p >= pr->points->size)
					break;
				ts = pr->points->times[p];
				/* fprintf (stdout, " %d times[%d] = %lf [%d, %d] += %d\n", 
					t, p, ts, tfirst, tlast, tstep); */
				p++;
				if (ts > t + tstep)
					break;
				scales_tot[i] += pr->points->points[p-1];
				scales_sam[i]++;
			}
			p--;

			if (scales_sam[i] > 0) {
				// scales_sam[i] = 1;
				size = scales->scales[i].size;
				if (agg == OFF)
					res = (double) scales_tot[i] / (double) scales_sam[i];
				else if (agg == ON)
					res = (double) scales_tot[i];
				scales->scales[i].points = (double*) util_MemRealloc (
					scales->scales[i].points, (size+1) * sizeof (double));
				scales->scales[i].points[size] = res;
				scales->scales[i].size++;
				scales_tot[i] = 0.0;
				scales_sam[i] = 0;
			}
		}
	}
	pr->scales_on = ON;
	long int mem = proc_SizeOfProcess (pr)
				+ sizeof (scales_sam) * conf->no +
				+ sizeof (scales_tot) * conf->no +
				+ sizeof (scales_pos) * conf->no;
	util_MemWr (mem);
	util_MemFree (scales_tot);
	util_MemFree (scales_sam);
	util_MemFree (scales_pos);
}

void proc_CreateScalesTime (
	proc_Process 		*pr,
	proc_ScalesConfig 	*conf)
{
	if (pr == NULL)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScalesTime");
	if (proc_CheckPoints (pr->points) == ERR ||
		pr->scales_on == ON || 
		pr->points->times == NULL)
		io_PrintErr (ERR, "invalid parameters for process %s in"
			" proc_CreateScalesTime", pr->name);
	if (proc_CheckScalesConfig (conf) == ERR)
		io_PrintErr (ERR, "invalid scale type for process %s in"
			" proc_CreateScalesTime", pr->name);

	proc_PrintHeader (pr);
	proc_PrintHeaderScales (conf);

	clock_t c;
	util_TimeIt (&c);

	proc_CreateScalesTimeAux (pr, conf, OFF);

	util_TimeWr (&c);
}


void proc_CreateScales (
	proc_Process *pr,
	proc_ScalesConfig *conf)
{
	if (pr == NULL)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScales");
	if (proc_CheckPoints (pr->points) == ERR ||
		pr->scales_on == ON)
		io_PrintErr (ERR, "invalid parameters for process %s in"
			" proc_CreateScales", pr->name);
	if (proc_CheckScalesConfig (conf) == ERR)
		io_PrintErr (ERR, "invalid scale type for process %s in"
			" proc_CreateScales", pr->name);

	proc_PrintHeader (pr);
	proc_PrintHeaderScales (conf);

	clock_t c;
	util_TimeIt (&c);

	int i, p;
	long long int size;
	double res;
	int *scales_sam = NULL;
	double *scales_tot = NULL;
	proc_Scales *scales = NULL;

	pr->scales = (proc_Scales*) util_MemMalloc (sizeof (proc_Scales));
	scales = pr->scales;
	scales->conf = conf;
	scales->scales = (proc_Points*) util_MemMalloc (conf->no * 
		sizeof (proc_Points));
	
	scales_sam = (int*)    util_MemCalloc (conf->no, sizeof (int));
	scales_tot = (double*) util_MemCalloc (conf->no, sizeof (double));
	
	for (i=0; i<conf->no; i++) {
		size = (int) ceil ((double) pr->points->size / conf->scales[i]);
		scales->scales[i].size = 0;
		scales->scales[i].scale = conf->scales[i];
		scales->scales[i].typeOfSignal = TestH_fGn;
		scales->scales[i].times = NULL; 

		scales->scales[i].points = (double*) util_MemCalloc (size, 
			sizeof (double));
	}

	for (p=0; p<pr->points->size; p++) {
		for (i=0; i<conf->no; i++) {
			scales_tot[i] += pr->points->points[p];
			scales_sam[i]++;
			if (scales_sam[i] == conf->scales[i] || 
				p == pr->points->size - 1) {
				
				size = scales->scales[i].size;
				res = (double) scales_tot[i] / (double) scales_sam[i];
				scales->scales[i].points[size] = res;
				scales->scales[i].size++;
				scales_tot[i] = 0.0;
				scales_sam[i] = 0;
			}
		}
	}
	pr->scales_on = ON;

	long int mem = proc_SizeOfProcess (pr)
				+ sizeof (scales_sam) * conf->no +
				+ sizeof (scales_tot) * conf->no;
	util_MemWr (mem);
	util_MemFree (scales_tot);
	util_MemFree (scales_sam);
	util_TimeWr (&c);
}

void proc_Aggregate (
	proc_Process *pr,
	int 		 m)
{
	if (pr == NULL || m <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScales");
	if (proc_CheckPoints (pr->points) == ERR ||
		pr->scales_on == ON)
		io_PrintErr (ERR, "invalid parameters for process %s in"
			" proc_CreateScales", pr->name);

	clock_t c;
	util_TimeIt (&c);

	int scale = pr->points->scale;
	tosig sig = pr->points->typeOfSignal;
	verb bak = TestHVerbosity;
	
	TestHVerbosity = TestH_NONE;
	proc_ScalesConfig *conf = proc_CreateScalesConfig (TestH_INC, m, m, 1);
	TestHVerbosity = bak;
	proc_CreateScalesTimeAux (pr, conf, ON);

	proc_DeletePoints (pr->points);
	pr->points = proc_CreatePoints (pr->scales->scales[0].points, NULL, 
		pr->scales->scales[0].size, scale, sig);
	proc_DeleteScales (pr);

	util_TimeWr (&c);
}

proc_Process* proc_Sum (
	proc_Process *pr,
	...) 
{



	return NULL;
}

void proc_Normalize (
	proc_Process *pr)
{
	if (pr == NULL || proc_CheckPoints (pr->points) == ERR)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" proc_NormalizeProcess");
	proc_Points *pts = pr->points;
	stat_Mean (pts);
	stat_StdDeviation (pts);
	int p;
	for (p=0; p<pts->size; p++) {
		pts->points[p] = (pts->points[p] - pts->mean) 
			/ pts->stdDev;
	}
	pr->points->typeOfSignal = TestH_fGn;
}

void proc_FractionalGaussianNoise (
	proc_Process *pr)
{
	if (pr == NULL || proc_CheckPoints (pr->points) == ERR || 
		pr->points->typeOfSignal != TestH_fBm)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" proc_FractionalGaussianNoise");
	proc_Points *pts = pr->points;
	int p;
	for (p=pts->size-1; p>0; p--)
		pts->points[p] -= pts->points[p - 1];
	pts->typeOfSignal = TestH_fGn;
}

void proc_FractionalBrownianMotion (
	proc_Process *pr)
{
	if (pr == NULL || proc_CheckPoints (pr->points) == ERR || 
		pr->points->typeOfSignal != TestH_fGn)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" proc_FractionalBrownianMotion");
	proc_Points *pts = pr->points;
	int p;
	for (p=1; p<pts->size; p++)
		pts->points[p] += pts->points[p - 1];
	pts->typeOfSignal = TestH_fBm;
}

/* void proc_CopyProcess (
	proc_Process *src,
	proc_Process *dst)
{
	if (src == NULL || dst != NULL ||
		proc_CheckPoints (src->points) == ERR)
		io_PrintErr (ERR, "invalid proc_Process structure(s) in"
			" proc_CopyProcess");

}

void proc_CopyPoints (
	proc_Points *src,
	proc_Points *dst)
{
	if (src == NULL || dst != NULL ||
		proc_CheckPoints (src) == ERR)
		io_PrintErr (ERR, "invalid proc_Points structure(s) in"
			" proc_CopyPoints");
	dst = create
} */

void proc_DeleteProcess ( 
	proc_Process *pr)
{
	if (pr != NULL) {
		proc_DeletePoints (pr->points);
		proc_DeleteScales (pr);
		pr = NULL;
	}
}

void proc_DeleteScales ( 
	proc_Process *pr)
{
	if (pr != NULL &&
		pr->scales_on == ON && pr->scales != NULL) {
		if (pr->scales->conf != NULL) {
			int s;
			for (s=0; s<pr->scales->conf->no; s++)
				proc_DeletePoints (&pr->scales->scales[s]);
			// util_MemFree (pr->scales->conf->scales);
		}
		util_MemFree (pr->scales->scales);
		util_MemFree (pr->scales);
		pr->scales_on = OFF;
		pr->scales = NULL;
	}
}

void proc_DeletePoints (
	proc_Points *pt)
{
	if (pt != NULL) {
		if (pt->points != NULL)
			util_MemFree (pt->points);
		if (pt->times != NULL)
			util_MemFree (pt->times);
		pt = NULL;
	}
}

int proc_CheckProcess ( 
	proc_Process *pr)
{
	if (proc_CheckPoints (pr->points) == ERR)
		return ERR;
	if (proc_CheckScales (pr->scales) == ERR)
		return ERR;
	return OK;
}

int proc_CheckScalesConfig (
	proc_ScalesConfig *conf)
{
	if (conf == NULL)
		return ERR;
	if (conf->typeOfScale != TestH_INC && conf->typeOfScale != TestH_POW &&
		conf->typeOfScale != TestH_RAND)
		return ERR;
	if (conf->min <= 0 || conf->max <= 0 || conf->min > conf->max || 
		conf->mov <= 0)
		return ERR;
	int i;
	for (i=0; i<conf->no; i++)
		if (conf->scales[i] < 1)
			return ERR;
	if (conf->typeOfScale == TestH_INC &&
		(conf->max - conf->min) / conf->mov + 1 < 1)
		return ERR;
	if (conf->typeOfScale == TestH_POW &&
		(conf->max - conf->min + 1 < 2 || conf->min < 1 || conf->mov < 2))
		return ERR;

	/* else if (typeOfScale == TestH_RAND && 
	   (scale_min < 2 || scale_max < 3 || scale_min <= scale_max))
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScales"); */

	return OK;
}

int proc_CheckScales ( 
	proc_Scales *sl)
{
	if (sl == NULL || sl->scales == NULL)
		return ERR;
	if (proc_CheckScalesConfig (sl->conf) == ERR)
		return ERR;
	int s;
	for (s=0; s<sl->conf->no; s++)
		if (proc_CheckPoints (&sl->scales[s]) == ERR)
			return ERR;
	return OK;
}

int proc_CheckPoints (
	proc_Points *pt)
{
	if (pt == NULL)
		return ERR;
	if (pt->points == NULL ||
		pt->size <= 0 || pt->scale <= 0 ||
		(pt->typeOfSignal != TestH_fGn && pt->typeOfSignal != TestH_fBm))
		return ERR;
	return OK;
}

void proc_PrintProcessStruct (
	proc_Process *pr)
{
	proc_PrintHeader (pr);
	if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF) {
		fprintf (stdout, "\n proc_Process ");
		TESTH_printVar (pr);
		if (pr == NULL)
			fprintf (stdout, " = NULL\n");
		else {
			fprintf (stdout, " {\n");
			fprintf (stdout, " \tname = \"%s\"\n", pr->name);
			if (pr->points == NULL) 
				fprintf (stdout, " \tproc_Points *points = NULL\n");
			else {
				fprintf (stdout, " \tproc_Points points {\n");
				if (pr->points->points == NULL)
					fprintf (stdout, " \t\tdouble *points = NULL\n");
				else
					fprintf (stdout, " \t\tdouble *points != NULL\n");
				if (pr->points->times == NULL)
					fprintf (stdout, " \t\tdouble *times = NULL\n");
				else
					fprintf (stdout, " \t\tdouble *times != NULL\n");
				if (pr->points->size <= 0)
					fprintf (stdout, " \t\tint size = undefined\n");
				else
					fprintf (stdout, " \t\tint size = %lld\n", 
						pr->points->size);
				if (pr->points->scale <= 0)
					fprintf (stdout, " \t\tint scale = undefined\n");
				else
					fprintf (stdout, " \t\tint scale = %d\n", 
						pr->points->scale);
				if (pr->points->typeOfSignal != TestH_fGn && 
					pr->points->typeOfSignal != TestH_fBm)
					fprintf (stdout, " \t\ttosig typeOfSignal = undefined\n");
				else if (pr->points->typeOfSignal == TestH_fGn)
					fprintf (stdout, " \t\ttosig typeOfSignal = TestH_fGn\n");
				else if (pr->points->typeOfSignal == TestH_fBm)
					fprintf (stdout, " \t\ttosig typeOfSignal = TestH_fBm\n");
				fprintf (stdout, " \t}\n");
			}
			if (pr->scales == NULL) 
				fprintf (stdout, " \tproc_Scales *scales = NULL\n");
			else {
				fprintf (stdout, " \tproc_Scales scales {\n");

				if (pr->scales->scales == NULL)
					fprintf (stdout, " \t\tproc_Points *scales = NULL\n");
				else
					fprintf (stdout, " \t\tproc_Points *scales != NULL\n");
				if (pr->scales->conf == NULL)
					fprintf (stdout, " \t\tproc_ScalesConfig *conf = NULL\n");
				else {
					fprintf (stdout, " \t\tproc_ScalesConfig *conf {\n");

					fprintf (stdout, "\t\t\ttos typeOfScale = %s\n", 
						proc_TypeOfScale (pr->scales->conf->typeOfScale));
					if (pr->scales->conf->scales == NULL)
						fprintf (stdout, " \t\t\tint *scales = NULL\n");
					else
						fprintf (stdout, " \t\t\tint *scales != NULL\n");
					if (pr->scales->conf->no <= 0)
						fprintf (stdout, " \t\t\tint no = undefined\n");
					else
						fprintf (stdout, " \t\t\tint no = %d\n", 
							pr->scales->conf->no);
					if (pr->scales->conf->min <= 0)
						fprintf (stdout, " \t\t\tint min = undefined\n");
					else
						fprintf (stdout, " \t\t\tint min = %d\n", 
							pr->scales->conf->min);
					if (pr->scales->conf->max <= 0)
						fprintf (stdout, " \t\t\tint max = undefined\n");
					else
						fprintf (stdout, " \t\t\tint max = %d\n", 
							pr->scales->conf->max);
					if (pr->scales->conf->mov <= 0)
						fprintf (stdout, " \t\t\tint mov = undefined\n");
					else
						fprintf (stdout, " \t\t\tint mov = %d\n", 
							pr->scales->conf->mov);
					fprintf (stdout, " \t\t}\n");
				}
				fprintf (stdout, " \t}\n");
			}
			fprintf (stdout, " }\n");
		}
		fprintf (stdout, "\n");
	}
}

void proc_PrintHeader (
	proc_Process *pr)
{
	if (TestHVerbosity > TestH_NONE &&
		TestHPrintPlain == OFF && TestHPrintHeader == ON) {
		if (pr == NULL || proc_CheckPoints (pr->points) == ERR)
			return;

		io_PrintSep ();        
		fprintf (stdout, "PROCESS: %s\n - %lld points",
			pr->name, pr->points->size);
		if (pr->points->typeOfSignal == TestH_fGn)
			fprintf (stdout, ", fGn\n");
		else if (pr->points->typeOfSignal == TestH_fBm)
			fprintf (stdout, ", fBm\n");
		else
			fprintf (stdout, "\n");
		if (pr->scales_on == ON)
			proc_PrintHeaderScales (pr->scales->conf);
			/* fprintf (stdout, " - %d scales, min: %d, max: %d, increment: %d\n",
				pr->scales->scale_no,
				pr->scales->scale_min, pr->scales->scale_max,
				pr->scales->mov); */
	}
}

void proc_PrintProcess (
	proc_Process *pr)
{
	proc_PrintHeader (pr);
	if (TestHVerbosity > TestH_NONE && pr->scales_on == ON) {
		int s;
		fprintf (stdout, "\n");
		if (TestHPrintPlain == OFF)
			fprintf (stdout, "\t scale\t\t points\n");
		for (s=0; s<pr->scales->conf->no; s++) {
			fprintf (stdout, " %12d\t%15lld\n",
				pr->scales->scales[s].scale, pr->scales->scales[s].size);
		}
		fprintf (stdout, "\n");
	}
}

void proc_PrintScales (
	proc_Process *pr)
{
	proc_PrintHeader (pr);
	if (TestHVerbosity > TestH_NONE) {
		proc_Scales *scales = pr->scales;
		if (proc_CheckScales (scales) == ERR)
			return;
		proc_ScalesConfig *conf = scales->conf;
		if (proc_CheckScalesConfig (conf) == ERR)
			return;

		int p, s;
		fprintf (stdout, "\n");
		if (TestHPrintPlain == OFF) {
			for (s=0; s<conf->no; s++)
				fprintf (stdout, "%12d    ", conf->scales[s]);
			fprintf (stdout, "\n");
			for (s=0; s<conf->no; s++)
				fprintf (stdout, "%12lld    ", scales->scales[s].size);
			fprintf (stdout, "\n");
		}
		for (p=0; p<scales->scales[0].size; p++) {
			fprintf (stdout, " ");
			for (s=0; s<scales->conf->no; s++)
				if (p < scales->scales[s].size)
					fprintf (stdout, "%15.8lf ", scales->scales[s].points[p]);
			fprintf (stdout, "\n");
		}
		fprintf (stdout, "\n");
	}
}

void proc_PrintPoints (
	proc_Points *pt)
{
	if (TestHVerbosity > TestH_NONE) {
		if (proc_CheckPoints (pt) == ERR)
			return;
		int p;
		fprintf (stdout, "\n");
		if (TestHPrintPlain == OFF) {
			fprintf (stdout, "\tindex");
			if (pt->times != NULL) {
				fprintf (stdout, "\t\t   times\t     points\n");
			} else
				fprintf (stdout, "\t\t points\n");
		}
		for (p=0; p<pt->size; p++) {
			fprintf (stdout, "  %10d", p+1);
			if (pt->times != NULL)
				fprintf (stdout, " %20.4lf", pt->times[p]);
			fprintf (stdout, " %20.6lf\n", pt->points[p]);

		}
		fprintf (stdout, "\n");
	}
}

void proc_PrintHeaderScales (
	proc_ScalesConfig *conf)
{
	if (TestHVerbosity > TestH_NONE && TestHPrintPlain == OFF) {
		fprintf (stdout, "TypeOfScale: %s\n", 
			proc_TypeOfScale (conf->typeOfScale));
		fprintf (stdout, " - %d scales, min: %d, max: %d", 
			conf->no, conf->min, conf->max);
		if (conf->typeOfScale == TestH_INC)
			fprintf (stdout, ", increment: %d\n", conf->mov);
		else if (conf->typeOfScale == TestH_RAND)
			fprintf (stdout, "\n");
		else if (conf->typeOfScale == TestH_POW)
			fprintf (stdout, ", exp base: %d\n", conf->mov);
		int i;
		fprintf (stdout, "\n [");
		for (i=0; i<conf->no; i++) {
			fprintf (stdout, " %d", conf->scales[i]);
			if (i != conf->no - 1)
				fprintf (stdout, ",");
		// else if (typeOfScale == TestH
		}
		fprintf (stdout, "]\n\n");
	}
}

static char* proc_TypeOfScale ( 
	tos scl)
{
	if (scl == TestH_INC)
		return "TestH_INC"; //  - incremental scaling";
	else if (scl == TestH_RAND)
		return "TestH_RAND"; // - random scaling";
	else if (scl == TestH_POW)
		return "TestH_POW"; // - power scaling";
	else
		return "undefined";
}


long int proc_SizeOfProcess (
	proc_Process *pr)
{
	long int sum = 0;
	if (pr != NULL) {
		if (pr->points != NULL) {
			sum += sizeof (pr->points->points) * 
				pr->points->size;
			if (pr->points->times != NULL)
				sum += sizeof (pr->points->times) * 
					pr->points->size;
		}
		if (pr->scales != NULL && pr->scales->conf != NULL) {
			int s;
			for (s=0; s<pr->scales->conf->no; s++)
				sum += sizeof (pr->scales->scales[s].points) * 
					pr->scales->scales[s].size;
		}
	}
	return sum;
}
