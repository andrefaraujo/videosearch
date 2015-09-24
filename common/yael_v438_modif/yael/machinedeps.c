/*
Copyright © INRIA 2009-2014.
Authors: Matthijs Douze & Herve Jegou
Contact: matthijs.douze@inria.fr  herve.jegou@inria.fr

This file is part of Yael.

    Yael is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Yael is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Yael.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>

#include <sys/time.h>


#include "machinedeps.h"


static int count_cpu_from_env() {
  int ncpu;
   
  if(!getenv("YAEL_COUNT_CPU")) return 0; 
  
  if(sscanf(getenv("YAEL_COUNT_CPU"), "%d", &ncpu) != 1 || ncpu <= 0) {
    fprintf(stderr, "could not parse YAEL_CPU_COUNT environment variable, using default\n"); 
    return 0; 
  } 
  return ncpu; 
}

#ifdef __linux__

#define __USE_GNU
#include <sched.h>

int count_cpu (void)
{
  int ncpu = count_cpu_from_env(); 
  if(ncpu) return ncpu; 

  cpu_set_t set;
  sched_getaffinity (0, sizeof (cpu_set_t), &set);
  int i, count = 0;
  for (i = 0; i < CPU_SETSIZE; i++)
    if (CPU_ISSET (i, &set))
      count++;
  return count;
}


#elif defined(__APPLE__)

#include <sys/types.h>
#include <sys/sysctl.h>


int count_cpu (void) {
  int ncpu = count_cpu_from_env(); 
  if(ncpu) return ncpu; 

  int count=-1;
  size_t count_size=sizeof(count);
  sysctlbyname("hw.ncpu",&count,&count_size,NULL,0);
  return count;
}

#else

int count_cpu() {
  return 1;
}


#endif

#ifndef __APPLE__

double log2(double x) {
  return log(x)/M_LN2;
}


#endif

#ifndef __linux__
void *memalign (size_t ignored, size_t nbytes)
{
  return malloc (nbytes);
}
#endif





double getmillisecs() 
{
  struct timeval tv;
  gettimeofday (&tv,NULL);
  return tv.tv_sec*1e3 +tv.tv_usec*1e-3;
}


/***********************************************************************
 *           Implementation of the threading part
 *
 * generic thread stuff */

#ifdef _OPENMP 

#include <omp.h>


#define GET_THREAD_NUM omp_get_thread_num()

#else


#define GET_THREAD_NUM 0

/* #pragma's will be ignored */

#endif


void compute_tasks (int n, int nt,
                    void (*task_fun) (void *arg, int tid, int i),
                    void *task_arg)
{
  int i;

#pragma omp parallel for schedule(dynamic) num_threads(nt)
  for(i = 0; i < n; i++) 
    (*task_fun)(task_arg, GET_THREAD_NUM, i);

}
