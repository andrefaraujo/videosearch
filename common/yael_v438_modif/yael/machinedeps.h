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

#ifndef MACHINEDEPS_H_INCLUDED
#define MACHINEDEPS_H_INCLUDED

#ifdef __APPLE__
#define HAVE_QSORT_R
#endif

#ifdef __linux__
#define HAVE_TLS
#else
#define __thread 
#endif

/*---------------------------------------------------------------------------*/
/*! @addtogroup machinedep
 *  @{
 */


/*! Return the number of cores. */
int count_cpu(void);

#ifndef __APPLE__

double log2(double x);

#endif

#ifdef __linux__
#include <malloc.h>
#else
#include <stdlib.h>

/*! allocate memory such that the pointer is aligned*/
void *memalign (size_t ignored, size_t nbytes);

#endif






/*! return a timestamp, which is useful to measure elapsed time */
double getmillisecs();

/*! exectutes a set of tasks in parallel using a thread pool 
 *
 * @param n            number of tasks to execute
 * @param nthread      number of threads that will run the tasks
 * @param task_fun     this callback will be called with 
 *              -  arg = task_arg
 *              -  tid = identifier of the thread in 0..nthread-1
 *              -  i = call number in 0..n-1
*/
void compute_tasks (int n, int nthread,
                    void (*task_fun) (void *arg, int tid, int i),
                    void *task_arg);




/*! @} */
#endif
