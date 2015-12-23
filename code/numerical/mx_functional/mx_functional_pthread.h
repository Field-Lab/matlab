/* mx_functional_pthread.h
 * Functional library for parallel Matlab MEX functions
 * Ver 0.1.0
 * Peter H. Li 2013 FreeBSD License
 *
 * See mx_functional for general documentation of the different functors.
 * The other classes here (messages and actors) are there to support the 
 * functors with multithreading capability based on pthread.
 * 
 */

#ifndef MX_FUNCTIONAL_PTHREAD
#define MX_FUNCTIONAL_PTHREAD

#include "mex.h"
#include "mx_functional.h"
#include <pthread.h>
#include <string.h>
#include <unistd.h>


/*
 * Generic column apply pthread functor template
 */
template <typename F, typename R = MxFunctionalNoOutput>
class ColApplyPthreadFunct
{
public:
  ColApplyPthreadFunct(F f, MWArrayContainer<R> out, int nthreads) : 
    _f(f), _out(out), _nthreads(nthreads) {}
  
  template <typename D>
  void *operator()(MWArrayContainer<D> a) {
    mwIndex ncols = a.ncols();
    
    for (mwIndex col = 0; col < ncols; col++) {
    }
    
    for (int col = 0; col < ncols; ++col);
  }
  
private:
  F _f;
  MWArrayContainer<R> _out;
  int _nthreads;
};


struct capt_args_t {
  int i;
  void *capt_funct;
  void *a;
};


template <typename F>
class ColApplyPthreadFunct<F,MxFunctionalNoOutput>
{
public:
  ColApplyPthreadFunct(F f, int nthreads) : 
    _f(f), _nthreads(nthreads) {}
  
  template <typename D>
  void operator()(MWArrayContainer<D> a) {
    mwIndex ncols = a.ncols();

    capt_args_t *argss = reinterpret_cast<capt_args_t *>(malloc(ncols * sizeof (capt_args_t)));
    memset(argss, 0, ncols * sizeof (capt_args_t));
    pthread_t *thread_ids = reinterpret_cast<pthread_t *>(malloc(ncols * sizeof (pthread_t)));
    int rc;
    for (mwIndex col = 0; col < ncols; ++col) {
      argss[col].i = col;
      argss[col].capt_funct = this;
      argss[col].a = &a;
      
      rc = pthread_create(&thread_ids[col], NULL, ColApplyPthreadFunct::run<D>, &argss[col]);
    }

    // Wait for threads to finish
    for (int col = 0; col < ncols; ++col) {
      rc = pthread_join(thread_ids[col], NULL);
    }
    
    // Cleanup
    free(argss);
    free(thread_ids);
  }
  
  template <typename D>
  void apply(MWArrayContainer<D> a, int col) {
    _f(a.get_col(col));
  }
  
private:
  template <typename D>
  static void *run(void *args_) {
    capt_args_t *args = reinterpret_cast<capt_args_t *>(args_);
    ColApplyPthreadFunct *capt_funct = reinterpret_cast<ColApplyPthreadFunct *>(args->capt_funct);
    MWArrayContainer<D> *a = reinterpret_cast<MWArrayContainer<D> *>(args->a);
    capt_funct->apply(*a, args->i);
  }
    
  F _f;
  int _nthreads;
};

#endif
