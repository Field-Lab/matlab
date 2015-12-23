/* mx_functional_parallel.h
 * Functional library for parallel Matlab MEX functions
 * Ver 0.3.0
 * Peter H. Li 2011 FreeBSD License
 *
 * See mx_functional for general documentation of the different functors.
 * The other classes here (messages and actors) are there to support the 
 * functors with multithreading capability based on the Theron Actors library.
 * See my blog post on getting Theron setup for Matlab:
 *   http://absurdlycertain.blogspot.com/2011/09/preamble-what-follows-is-guide.html
 */

#ifndef MX_FUNCTIONAL_PARALLEL
#define MX_FUNCTIONAL_PARALLEL

#include "mex.h"
#include "mx_functional.h"

#define THERON_USE_BOOST_THREADS 1
#include "Theron/Framework.h"
#include "Theron/Actor.h"
#include "Theron/Receiver.h"

template <class F, class D, class R = MxFunctionalNoOutput>
struct ColApplyMessage
{
  F f;
  MWArrayContainer<D> in;
  MWArrayContainer<R> out;
};

template <class F, class D>
struct ColApplyMessage<F,D,MxFunctionalNoOutput>
{
  F f;
  MWArrayContainer<D> in;
};


template <class F, class D, class R = MxFunctionalNoOutput>
class ColApplyActor : public Theron::Actor
{
public:
  ColApplyActor() {
    RegisterHandler(this, &ColApplyActor<F,D,R>::RequestHandler);
  }
  
private:
  void RequestHandler(const ColApplyMessage<F,D,R> &message, const Theron::Address from) {
    message.f(message.in, message.out);
    TailSend(0, from);
  }
};

template <class F, class D>
class ColApplyActor<F,D,MxFunctionalNoOutput> : public Theron::Actor
{
public:
  ColApplyActor() {
    RegisterHandler(this, &ColApplyActor<F,D>::RequestHandler);
  }
  
private:
  void RequestHandler(const ColApplyMessage<F,D> &message, const Theron::Address from) {
    message.f(message.in);
    TailSend(0, from);
  }
};


template <typename F, typename R = MxFunctionalNoOutput>
class ColApplyParallelFunct
{
public:
  ColApplyParallelFunct(F f, MWArrayContainer<R> out, int nthreads) : 
    _f(f), _out(out), _nthreads(nthreads) {}
  
  template <typename D>
  void *operator()(MWArrayContainer<D> a) {
    mwIndex ncols = a.ncols();
      
    Theron::Framework theron(_nthreads);
    Theron::Receiver receiver;
    
    Theron::ActorRef actors[ncols];
    for (mwIndex col = 0; col < ncols; col++) {
      ColApplyMessage<F,D,R> message = {_f, a.get_col(col), _out.get_col(col)};
      actors[col] = Theron::ActorRef(theron.CreateActor< ColApplyActor<F,D,R> >());
      actors[col].Push(message, receiver.GetAddress());
    }
    
    for (int col = 0; col < ncols; ++col) receiver.Wait();
  }
  
private:
  F _f;
  MWArrayContainer<R> _out;
  int _nthreads;
};


template <typename F>
class ColApplyParallelFunct<F,mxArray *>
{
public:
  ColApplyParallelFunct(F f, mxArray *out, int nthreads) : 
    _f(f), _out(out), _nthreads(nthreads) {}
  
  template <typename D>
  void *operator()(MWArrayContainer<D> a) {
    mwIndex ncols = a.ncols();
      
    Theron::Framework theron(_nthreads);
    Theron::Receiver receiver;
    
    // Assume cast output to MWArrayContainer<D> is desired
    MWArrayContainer<D> out(_out);
    
    Theron::ActorRef actors[ncols];
    for (mwIndex col = 0; col < ncols; ++col) {
      ColApplyMessage<F,D,D> message = {_f, a.get_col(col), out.get_col(col)};
      actors[col] = Theron::ActorRef(theron.CreateActor< ColApplyActor<F,D,D> >());
      actors[col].Push(message, receiver.GetAddress());
    }
    
    for (int col = 0; col < ncols; ++col) receiver.Wait();
  }
  
private:
  F _f;
  mxArray *_out;
  int _nthreads;
};


template <typename F>
class ColApplyParallelFunct<F,MxFunctionalNoOutput>
{
public:
  ColApplyParallelFunct(F f, int nthreads) : 
    _f(f), _nthreads(nthreads) {}
  
  template <typename D>
  void operator()(MWArrayContainer<D> a) {
    mwIndex ncols = a.ncols();
    
    Theron::Framework theron(_nthreads);
    Theron::Receiver receiver;

    Theron::ActorRef actors[ncols];
    for (mwIndex col = 0; col < ncols; ++col) {
      ColApplyMessage<F,D> message = {_f, a.get_col(col)};
      actors[col] = Theron::ActorRef(theron.CreateActor< ColApplyActor<F,D> >());
      actors[col].Push(message, receiver.GetAddress());
    }
    
    for (int col = 0; col < ncols; ++col) receiver.Wait();
  }
  
private:
  F _f;
  int _nthreads;
};






template <class F>
struct BlockApplyMessage
{
  F f;
  int blocknum;
  int nblocks;
};


template <class F>
class BlockApplyActor : public Theron::Actor
{
public:
  BlockApplyActor() {
    RegisterHandler(this, &BlockApplyActor<F>::RequestHandler);
  }
  
private:
  void RequestHandler(const BlockApplyMessage<F> &message, const Theron::Address from) {
    message.f(message.blocknum, message.nblocks);
    TailSend(0, from);
  }
};


template <typename F>
class BlockApplyParallelFunct
{
public:
  BlockApplyParallelFunct(F f) : _f(f) {}
  
  void operator()(int nblocks) {
    Theron::Framework theron(nblocks);
    Theron::Receiver receiver;

    Theron::ActorRef actors[nblocks];
    for (int blk = 0; blk < nblocks; ++blk) {
      BlockApplyMessage<F> message = {_f, blk, nblocks};
      actors[blk] = Theron::ActorRef(theron.CreateActor< BlockApplyActor<F> >());
      actors[blk].Push(message, receiver.GetAddress());
    }
    
    for (int blk = 0; blk < nblocks; ++blk) receiver.Wait();
  }
private:
  F _f;
};


#endif
