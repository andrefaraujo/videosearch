/*
Copyright © INRIA 2010-2011. 
Authors: Matthijs Douze & Herve Jegou 
Contact: matthijs.douze@inria.fr  herve.jegou@inria.fr

This software is a computer program whose purpose is to provide 
efficient tools for basic yet computationally demanding tasks, 
such as find k-nearest neighbors using exhaustive search 
and kmeans clustering. 

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
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


/*! trace all mallocs between two function calls. Intended to replace
 * struct mallinfo that does not seem to work. Implemented only for
 * Linux. Includes inefficient code that should not be relied on while
 * profiling. */

typedef struct {
  int n_alloc,n_free,n_realloc; /* nb of operations of each type */
  size_t delta_alloc;           /* total allocated minus deallocated 
                                   (can only trace deallocs that were allocated during the tracing) */ 
  size_t max_alloc;             /* max of delta_alloc during the tracing */
  int n_untracked_frees;        /* nb of frees of objects that were not allocated during the tracing 
                                   (cannot be taken into accout in delta_alloc) */
} malloc_stats_t;


void malloc_stats_begin (void);

malloc_stats_t malloc_stats_end (void);



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
