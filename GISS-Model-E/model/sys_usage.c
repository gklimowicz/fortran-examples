#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <sys/time.h>
#include <sys/resource.h>

#ifndef RUSAGE_SELF
#define RUSAGE_SELF 0
#endif


/*
void PrintUsage()
{
struct rusage buf;
static float u_time=0., s_time=0.;

getrusage( RUSAGE_SELF, &buf );

printf("#################### SYSTEM USAGE ####################\n");
printf("User Time: %f\n",   (float)buf.ru_utime.tv_sec + (float)buf.ru_utime.tv_usec*1.e-6 - u_time);
printf("System Time: %f\n", (float)buf.ru_stime.tv_sec + (float)buf.ru_stime.tv_usec*1.e-6 - s_time);

u_time = (float)buf.ru_utime.tv_sec + (float)buf.ru_utime.tv_usec*1.e-6;
s_time = (float)buf.ru_stime.tv_sec + (float)buf.ru_stime.tv_usec*1.e-6;

printf("Max Resident Memory: %d K\n", buf.ru_maxrss);
printf("Text Size: %d K\n", buf.ru_ixrss);
printf("Data Size: %d K\n", buf.ru_idrss);
printf("Stack Size: %d K\n", buf.ru_isrss);
printf("######################################################\n");

}
*/


void SysUsage(int var, int action)
{
struct rusage buf;
struct timeval buf1;
static float total[16], beg[16];
static float s_total[16], s_beg[16];
static double g_total[16], g_beg[16];
float u_time=0., s_time=0.;
static double g_time=0;
/* next line is for SunOS only */
int getrusage(int who, struct rusage *rusage);

getrusage( RUSAGE_SELF, &buf );
gettimeofday(&buf1, NULL);

u_time = (float)buf.ru_utime.tv_sec + (float)buf.ru_utime.tv_usec*1.e-6;
s_time = (float)buf.ru_stime.tv_sec + (float)buf.ru_stime.tv_usec*1.e-6;
/*g_time = time(NULL);*/
 g_time = (double)buf1.tv_sec + (double)buf1.tv_usec*1.e-6;

switch(action)
        {
        case 0: total[var] = 0.; s_total[var] = 0.; g_total[var] = 0.;
                break;

        case 1: beg[var] = u_time; s_beg[var] = s_time; g_beg[var] = g_time;
                break;

        case 2: total[var] += u_time - beg[var];
                s_total[var] += s_time - s_beg[var];
                g_total[var] += g_time - g_beg[var];
                break;

        case 3: printf("User Time for #%1d#   %f  System: %f Glob: %f\n",
                       var, total[var], s_total[var], (float)(g_total[var]));
                break;
        }

}

/*define fortran interface */
/*extern "C" {*/

void sysusage_(int *var, int *action)
{
  SysUsage( *var, *action);
}

/*}*/

