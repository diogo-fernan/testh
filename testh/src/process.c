
#include "process.h"
#include "io.h"
#include "utilities.h"
#include "regression.h"
#include "statistics.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static char* proc_TypeOfScale (tos typeOfScale);

proc_Process* proc_CreateProcess (
	const char 	*name,
	double 		*points,
	double 		*times,
	int 		size,
	tosig 		typeOfSignal)
{
	if (points == NULL || points == NULL ||
		size <= 0 || (typeOfSignal != TestH_fGn && typeOfSignal != TestH_fBm))
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateProcess");

	proc_Process *p = (proc_Process*) util_MemMalloc (sizeof (proc_Process));
	p->points = proc_CreatePoints (points, times, size, 1, typeOfSignal);
	p->scales = NULL;
	p->scales_on = OFF;
	if (name != NULL) {
		int s = strlen(name) + 1;
		p->name = (char*) util_MemMalloc (s * sizeof (char));
		memset (p->name, '\0', s);
		memcpy (p->name, name, s);
	} else {
		p->name = (char*) util_MemMalloc (11 * sizeof (char));
		memset (p->name, '\0', 11);
		p->name = "empty name";
	}

	return p;
}

proc_Points* proc_CreatePoints (
	double 	*points,
	double 	*times,
	int 	size,
	int 	scale,
	tosig 	typeOfSignal)
{
	if (points == NULL || size < 0 || scale < 0)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreatePoints");

	int p;
	proc_Points *pts = (proc_Points*) util_MemMalloc (sizeof (proc_Points));
	pts->size = size;
	pts->scale = scale;
	pts->typeOfSignal = typeOfSignal;
	pts->points = (double*) util_MemMalloc (size * sizeof (double));
	
	for (p=0; p<size; p++)
		pts->points[p] = points[p];
	
	if (times != NULL) {
		int t;
		pts->times = (double*) util_MemMalloc (size * sizeof (double));
		for (t=0; t<size; t++)
			pts->times[t] = times[t];
	} else
		pts->times = NULL;
		
	return pts;
}

// ... --> for manual scaling
proc_ScalesConfig* proc_CreateScalesConfig (
	tos typeOfScale,
	int min,
	int max,
	int mov,
	...)
{
	if (typeOfScale != TestH_INC && 
		typeOfScale != TestH_RAND && 
		typeOfScale != TestH_POW)
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
	conf->typeOfScale = typeOfScale;
	if (typeOfScale == TestH_INC)
		conf->no = (max - min) / mov + 1;
	// else if (typeOfScale == TestH_RAND)
	else // if (typeOfScale == TestH_POW) 
		conf->no = max - min + 1;

	io_PrintSep ();        
	// proc_PrintHeaderScales (conf);
	if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF) {
		fprintf (stdout, " - %d scales, min: %d, max: %d",
			conf->no, min, max);
		if (typeOfScale == TestH_INC)
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
		if (typeOfScale == TestH_POW)
			s = pow (mov, min + i);

		conf->scales[i] = s;
		if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF) {
			fprintf (stdout, " %lld", s);
			if (i != conf->no - 1)
				fprintf (stdout, ",");
		}

		if (typeOfScale == TestH_INC)
			s += mov;
		// else if (typeOfScale == TestH_RAND)
		// else if (typeOfScale == TestH_MAN)
	}
	if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF)
		fprintf (stdout, "]\n");

	if (proc_CheckScalesConfig (conf) == ERR)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScalesConfig");
	return conf;
}

void proc_CreateScalesTime (proc, scale_min, scale_max, scale_multiplier)
	proc_Process *proc;
	double scale_min;
	double scale_max;
	int scale_multiplier;
{
	/* if (proc_CheckPoints (proc->points) == ERR ||
		proc->scales_on == ON ||
		scale_min <= 0.0 || scale_max <= 0.0 || scale_multiplier <= 0)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScalesTime");

	int i, p, scale_no, size, flag;
	int *scales_sam, *scales_tot;
	double s, res;
	proc_Scales *scales;

	scale_no = 0;
	s = scale_min / scale_multiplier;
	while ((s *= scale_multiplier) <= scale_max)
		scale_no++;
		
	proc->scales = (proc_Scales*) util_MemMalloc (sizeof (proc_Scales));
	scales = proc->scales;
	scales->scale_no  = scale_no;
	scales->scale_min = scale_min;
	scales->scale_max = scale_max;
	scales->mov = scale_multiplier;
	scales->scales = (proc_Points*) util_MemMalloc (scale_no * 
		sizeof (proc_Points));
	
	scales_sam = (int*) util_MemCalloc (scale_no, sizeof (int));
	scales_tot = (int*) util_MemCalloc (scale_no, sizeof (int));

	i = 0;
	s = scale_min / scale_multiplier;
	while ((s *= scale_multiplier) <= scale_max) {
		size = 0;
		for (p=0; p<proc->points->size; p++) {
			if (proc->points->times[p] >= (double) (size + 1) * s ||
				p == proc->points->size - 1) {
				size++;
			}
		}
		scales->scales[i].size = 0;
		scales->scales[i].scale = s;
		scales->scales[i++].points = (double*) util_MemCalloc (size, 
			sizeof (double));
	}
	
	for (p=0; p<proc->points->size; p++) {
		i = 0;
		s = scale_min / scale_multiplier;
		while ((s *= scale_multiplier) <= scale_max) {
			flag = 0;
			if (proc->points->times[p+1] >=
					(double) (scales->scales[i].size + 1) * s) {
				flag = 1;
			} else if (p == proc->points->size - 1) {
				scales_tot[i] += proc->points->points[p];
				scales_sam[i]++;
				flag = 2;
			}
			
			if (flag) {
				if (flag == 1 && scales_tot[i] == 0) {
					scales_tot[i] += proc->points->points[p];
					scales_sam[i]++;
				}
				size = scales->scales[i].size;
				res = (double) scales_tot[i] / scales_sam[i];
				scales->scales[i].points[size] = res;
				scales->scales[i].size++;
				scales_tot[i] = 0;
				scales_sam[i] = 0;
				
			} else {
				scales_tot[i] += proc->points->points[p];
				scales_sam[i]++;
			}
			i++;
		}
	}
	
	proc->scales_on = ON;

	util_MemFree (scales_tot);
	util_MemFree (scales_sam); */
}


void proc_CreateScales (
	proc_Process *proc,
	proc_ScalesConfig *conf)
{
	if (proc == NULL)
		io_PrintErr (ERR, "invalid parameters in"
			" proc_CreateScales");
	if (proc_CheckPoints (proc->points) == ERR ||
		proc->scales_on == ON)
		io_PrintErr (ERR, "invalid parameters for process %s in"
			" proc_CreateScales", proc->name);
	if (proc_CheckScalesConfig (conf) == ERR)
		io_PrintErr (ERR, "invalid scale type for process %s in"
			" proc_CreateScales", proc->name);

	proc_PrintHeader (proc);
	proc_PrintHeaderScales (conf);

	clock_t c;
	util_TimeIt (&c);

	int i, p;
	long long int size;
	double res;
	int *scales_sam = NULL;
	double *scales_tot = NULL;
	proc_Scales *scales = NULL;

	proc->scales = (proc_Scales*) util_MemMalloc (sizeof (proc_Scales));
	scales = proc->scales;
	scales->conf = conf;
	scales->scales = (proc_Points*) util_MemMalloc (conf->no * 
		sizeof (proc_Points));
	
	scales_sam = (int*)    util_MemCalloc (conf->no, sizeof (int));
	scales_tot = (double*) util_MemCalloc (conf->no, sizeof (double));
	
	for (i=0; i<conf->no; i++) {
		size = (int) ceil ((double) proc->points->size / conf->scales[i]);
		scales->scales[i].size = 0;
		scales->scales[i].scale = conf->scales[i];
		scales->scales[i].typeOfSignal = TestH_fGn;
		scales->scales[i].times = NULL; 

		scales->scales[i].points = (double*) util_MemCalloc (size, 
			sizeof (double));
	}

	for (p=0; p<proc->points->size; p++) {
		for (i=0; i<conf->no; i++) {
			scales_tot[i] += proc->points->points[p];
			scales_sam[i]++;
			if (scales_sam[i] == conf->scales[i] || 
				p == proc->points->size - 1) {
				
				size = scales->scales[i].size;
				res = (double) scales_tot[i] / (double) scales_sam[i];
				scales->scales[i].points[size] = res;
				scales->scales[i].size++;
				scales_tot[i] = 0.0;
				scales_sam[i] = 0;
			}
		}
	}
	proc->scales_on = ON;

	long int mem = proc_SizeOfProcess (proc)
				+ sizeof (scales_sam) * conf->no +
				+ sizeof (scales_tot) * conf->no;
	util_MemWr (mem);
	util_MemFree (scales_tot);
	util_MemFree (scales_sam);
	util_TimeWr (&c);
}


void proc_Normalize (proc)
	proc_Process *proc;
{
	if (proc == NULL || proc_CheckPoints (proc->points) == ERR)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" proc_NormalizeProcess");
	proc_Points *pts = proc->points;
	stat_Mean (pts);
	stat_StdDeviation (pts);
	int p;
	for (p=0; p<pts->size; p++) {
		pts->points[p] = (pts->points[p] - pts->mean) 
			/ pts->stdDev;
	}
	proc->points->typeOfSignal = TestH_fGn;
}

void proc_FractionalGaussianNoise (proc)
	proc_Process *proc;
{
	if (proc == NULL || proc_CheckPoints (proc->points) == ERR || 
		proc->points->typeOfSignal != TestH_fBm)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" proc_FractionalGaussianNoise");
	proc_Points *pts = proc->points;
	int p;
	for (p=pts->size-1; p>0; p--)
		pts->points[p] -= pts->points[p - 1];
	pts->typeOfSignal = TestH_fGn;
}

void proc_FractionalBrownianMotion (proc)
	proc_Process *proc;
{
	if (proc == NULL || proc_CheckPoints (proc->points) == ERR || 
		proc->points->typeOfSignal != TestH_fGn)
		io_PrintErr (ERR, "invalid proc_Process structure in"
			" proc_FractionalBrownianMotion");
	proc_Points *pts = proc->points;
	int p;
	for (p=1; p<pts->size; p++)
		pts->points[p] += pts->points[p - 1];
	pts->typeOfSignal = TestH_fBm;
}


void proc_DeleteProcess (proc) 
	proc_Process *proc;
{
	if (proc != NULL) {
		proc_DeletePoints (proc->points);
		proc_DeleteScales (proc);
		proc = NULL;
	}
}

void proc_DeleteScales (proc) 
	proc_Process *proc;
{
	if (proc != NULL &&
		proc->scales_on == ON && proc->scales != NULL) {
		if (proc->scales->conf != NULL) {
			int s;
			for (s=0; s<proc->scales->conf->no; s++)
				proc_DeletePoints (&proc->scales->scales[s]);
			// util_MemFree (proc->scales->conf->scales);
		}
		util_MemFree (proc->scales->scales);
		util_MemFree (proc->scales);
		proc->scales_on = OFF;
		proc->scales = NULL;
	}
}

void proc_DeletePoints (pts)
	proc_Points *pts;
{
	if (pts != NULL) {
		if (pts->points != NULL)
			util_MemFree (pts->points);
		if (pts->times != NULL)
			util_MemFree (pts->times);
		pts = NULL;
	}
}

int proc_CheckProcess (proc) 
	proc_Process *proc;
{
	if (proc_CheckPoints (proc->points) == ERR)
		return ERR;
	if (proc_CheckScales (proc->scales) == ERR)
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
	if (conf->min <= 0 || conf->max <= 0 || conf->min >= conf->max || 
		conf->mov <= 0)
		return ERR;
	int i;
	for (i=0; i<conf->no; i++)
		if (conf->scales[i] < 2)
			return ERR;
	if (conf->typeOfScale == TestH_INC &&
		(conf->max - conf->min) / conf->mov + 1 < 2)
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

int proc_CheckScales (scales) 
	proc_Scales *scales;
{
	if (scales == NULL || scales->scales == NULL)
		return ERR;
	if (proc_CheckScalesConfig (scales->conf) == ERR)
		return ERR;
	int s;
	for (s=0; s<scales->conf->no; s++)
		if (proc_CheckPoints (&scales->scales[s]) == ERR)
			return ERR;
	return OK;
}

int proc_CheckPoints (pts)
	proc_Points *pts;
{
	if (pts == NULL)
		return ERR;
	if (pts->points == NULL ||
		pts->size <= 0 || pts->scale <= 0 ||
		(pts->typeOfSignal != TestH_fGn && pts->typeOfSignal != TestH_fBm))
		return ERR;
	return OK;
}

void proc_PrintProcessStruct (proc)
	proc_Process *proc;
{
	proc_PrintHeader (proc);
	if (TestHVerbosity > TestH_LOW && TestHPrintPlain == OFF) {
		fprintf (stdout, "\n proc_Process ");
		TESTH_printVar (proc);
		if (proc == NULL)
			fprintf (stdout, " = NULL\n");
		else {
			fprintf (stdout, " {\n");
			fprintf (stdout, " \tname = \"%s\"\n", proc->name);
			if (proc->points == NULL) 
				fprintf (stdout, " \tproc_Points *points = NULL\n");
			else {
				fprintf (stdout, " \tproc_Points points {\n");
				if (proc->points->points == NULL)
					fprintf (stdout, " \t\tdouble *points = NULL\n");
				else
					fprintf (stdout, " \t\tdouble *points != NULL\n");
				if (proc->points->times == NULL)
					fprintf (stdout, " \t\tdouble *times = NULL\n");
				else
					fprintf (stdout, " \t\tdouble *times != NULL\n");
				if (proc->points->size <= 0)
					fprintf (stdout, " \t\tint size = undefined\n");
				else
					fprintf (stdout, " \t\tint size = %lld\n", 
						proc->points->size);
				if (proc->points->scale <= 0)
					fprintf (stdout, " \t\tint scale = undefined\n");
				else
					fprintf (stdout, " \t\tint scale = %d\n", 
						proc->points->scale);
				if (proc->points->typeOfSignal != TestH_fGn && 
					proc->points->typeOfSignal != TestH_fBm)
					fprintf (stdout, " \t\ttosig typeOfSignal = undefined\n");
				else if (proc->points->typeOfSignal == TestH_fGn)
					fprintf (stdout, " \t\ttosig typeOfSignal = TestH_fGn\n");
				else if (proc->points->typeOfSignal == TestH_fBm)
					fprintf (stdout, " \t\ttosig typeOfSignal = TestH_fBm\n");
				fprintf (stdout, " \t}\n");
			}
			if (proc->scales == NULL) 
				fprintf (stdout, " \tproc_Scales *scales = NULL\n");
			else {
				fprintf (stdout, " \tproc_Scales scales {\n");

				if (proc->scales->scales == NULL)
					fprintf (stdout, " \t\tproc_Points *scales = NULL\n");
				else
					fprintf (stdout, " \t\tproc_Points *scales != NULL\n");
				if (proc->scales->conf == NULL)
					fprintf (stdout, " \t\tproc_ScalesConfig *conf = NULL\n");
				else {
					fprintf (stdout, " \t\tproc_ScalesConfig *conf {\n");

					fprintf (stdout, "\t\t\ttos typeOfScale = %s\n", 
						proc_TypeOfScale (proc->scales->conf->typeOfScale));
					if (proc->scales->conf->scales == NULL)
						fprintf (stdout, " \t\t\tint *scales = NULL\n");
					else
						fprintf (stdout, " \t\t\tint *scales != NULL\n");
					if (proc->scales->conf->no <= 0)
						fprintf (stdout, " \t\t\tint no = undefined\n");
					else
						fprintf (stdout, " \t\t\tint no = %d\n", 
							proc->scales->conf->no);
					if (proc->scales->conf->min <= 0)
						fprintf (stdout, " \t\t\tint min = undefined\n");
					else
						fprintf (stdout, " \t\t\tint min = %d\n", 
							proc->scales->conf->min);
					if (proc->scales->conf->max <= 0)
						fprintf (stdout, " \t\t\tint max = undefined\n");
					else
						fprintf (stdout, " \t\t\tint max = %d\n", 
							proc->scales->conf->max);
					if (proc->scales->conf->mov <= 0)
						fprintf (stdout, " \t\t\tint mov = undefined\n");
					else
						fprintf (stdout, " \t\t\tint mov = %d\n", 
							proc->scales->conf->mov);
					fprintf (stdout, " \t\t}\n");
				}
				fprintf (stdout, " \t}\n");
			}
			fprintf (stdout, " }\n");
		}
		fprintf (stdout, "\n");
	}
}

void proc_PrintHeader (proc)
	proc_Process *proc;
{
	if (TestHVerbosity > TestH_NONE &&
		TestHPrintPlain == OFF && TestHPrintHeader == ON) {
		if (proc == NULL || proc_CheckPoints (proc->points) == ERR)
			return;

		io_PrintSep ();        
		fprintf (stdout, "PROCESS: %s\n - %lld points",
			proc->name, proc->points->size);
		if (proc->points->typeOfSignal == TestH_fGn)
			fprintf (stdout, ", fGn\n");
		else if (proc->points->typeOfSignal == TestH_fBm)
			fprintf (stdout, ", fBm\n");
		else
			fprintf (stdout, "\n");
		if (proc->scales_on == ON)
			proc_PrintHeaderScales (proc->scales->conf);
			/* fprintf (stdout, " - %d scales, min: %d, max: %d, increment: %d\n",
				proc->scales->scale_no,
				proc->scales->scale_min, proc->scales->scale_max,
				proc->scales->mov); */
	}
}

void proc_PrintProcess (proc)
	proc_Process *proc;
{
	proc_PrintHeader (proc);
	if (TestHVerbosity > TestH_NONE && proc->scales_on == ON) {
		int s;
		fprintf (stdout, "\n");
		if (TestHPrintPlain == OFF)
			fprintf (stdout, "\tscale\t\t  points\n");
		for (s=0; s<proc->scales->conf->no; s++) {
			fprintf (stdout, " %12d\t%15lld\n",
				proc->scales->scales[s].scale, proc->scales->scales[s].size);
		}
		fprintf (stdout, "\n");
	}
}

void proc_PrintScales (proc)
	proc_Process *proc;
{
	proc_PrintHeader (proc);
	if (TestHVerbosity > TestH_NONE) {
		proc_Scales *scales = proc->scales;
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

void proc_PrintPoints (pts)
	proc_Points *pts;
{
	if (TestHVerbosity > TestH_NONE) {
		if (proc_CheckPoints (pts) == ERR)
			return;
		int p;
		fprintf (stdout, "\n");
		if (TestHPrintPlain == OFF)
			fprintf (stdout, "\tindex\t\t\t   points\n");
		for (p=0; p<pts->size; p++)
			fprintf (stdout, "  %10d\t\t%20.10lf\n", p+1, pts->points[p]);
		fprintf (stdout, "\n");
	}
}

void proc_PrintHeaderScales (proc_ScalesConfig *conf)
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
			// else if (typeOfScale == TestH_RAND)
		}
		fprintf (stdout, "]\n");
	}
}

static char* proc_TypeOfScale ( 
	tos typeOfScale)
{
	if (typeOfScale == TestH_INC)
		return "TestH_INC"; //  - incremental scaling";
	else if (typeOfScale == TestH_RAND)
		return "TestH_RAND"; // - random scaling";
	else if (typeOfScale == TestH_POW)
		return "TestH_POW"; // - power scaling";
	else
		return "undefined";
}


long int proc_SizeOfProcess (
	proc_Process *proc)
{
	long int sum = 0;
	if (proc != NULL) {
		if (proc->points != NULL) {
			sum += sizeof (proc->points->points) * 
				proc->points->size;
			if (proc->points->times != NULL)
				sum += sizeof (proc->points->times) * 
					proc->points->size;
		}
		if (proc->scales != NULL && proc->scales->conf != NULL) {
			int s;
			for (s=0; s<proc->scales->conf->no; s++)
				sum += sizeof (proc->scales->scales[s].points) * 
					proc->scales->scales[s].size;
		}
	}
	return sum;
}
