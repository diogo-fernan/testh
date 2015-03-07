#ifndef __TESTH_PROCESS__
#define __TESTH_PROCESS__

enum TestHTypeOfSignal { 
	TestH_fGn 	= 0x1,
	TestH_fBm 	= 0x2
};

enum TestHTypeOfScale {
	TestH_INC	= 0x1,
	TestH_RAND	= 0x1,
	TestH_POW	= 0x2
};

typedef enum TestHTypeOfSignal tosig;
typedef enum TestHTypeOfScale  tos;


typedef struct {
	tos typeOfScale;

	unsigned int min;
	unsigned int max;
	unsigned int no;
	unsigned int mov;

	int *scales;
} proc_ScalesConfig;



typedef struct {
	double			*points;
	double			*times;
	unsigned int	scale; // m
	unsigned long long int	size;

	tosig typeOfSignal;

	double mean;
	double stdDev;

	double amt;
	unsigned int moment;
} proc_Points;

typedef struct {
	proc_Points *scales;
	proc_ScalesConfig *conf;
} proc_Scales;

typedef struct {
	proc_Points *points;
	proc_Scales *scales; // malloc (sizeof (proc_Scales))
	int scales_on;
	char *name;
} proc_Process;



proc_Process* proc_CreateProcess (const char *name, double *points,
	double *times, int size, tosig typeOfSignal);
proc_Points* proc_CreatePoints (double *points, double *times,
	int size, int scale, tosig typeOfSignal);

proc_ScalesConfig* proc_CreateScalesConfig (tos typeOfScale, int min, int max, 
	int mov, ...);
void proc_CreateScales (proc_Process *proc, proc_ScalesConfig *conf);
void proc_CreateScalesTime (proc_Process *proc, double scale_min,
	double scale_max, int scale_multiplier);


void proc_Normalize (proc_Process *proc);
void proc_FractionalGaussianNoise (proc_Process *proc);
void proc_FractionalBrownianMotion (proc_Process *proc);

void proc_PrintProcessStruct (proc_Process *proc);

void proc_PrintHeader (proc_Process *proc);
void proc_PrintProcess (proc_Process *proc);
void proc_PrintScales (proc_Process *proc);
void proc_PrintPoints (proc_Points *points);

void proc_DeleteProcess (proc_Process *proc);
void proc_DeleteScales (proc_Process *proc);
void proc_DeletePoints (proc_Points *pts);

int proc_CheckPoints (proc_Points *points);
int proc_CheckScalesConfig (proc_ScalesConfig *conf);
int proc_CheckScales (proc_Scales *scales);
int proc_CheckProcess (proc_Process *proc);

void proc_PrintHeaderScales (proc_ScalesConfig *config);

long int proc_SizeOfProcess (proc_Process *proc);


#endif