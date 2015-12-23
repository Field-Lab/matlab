/* boxspline2_parallel.cpp
 * Version 0.8
 *
 * Parallelized version of boxspline2.  Depends on Theron Actors library
 * for cleaner multithreading.  See http://absurdlycertain.blogspot.com/2011/09/preamble-what-follows-is-guide.html
 * for a guide to compiling Mex functions with Theron.
 *
 * 2011-09 Peter H. Li
 * FreeBSD License
 */

// Compile hints
//
// See http://absurdlycertain.blogspot.com/2011/09/preamble-what-follows-is-guide.html
// for more information.
//
// Mac with macports:
// mex -DTHERON_USE_BOOST_THREADS -I/Users/peterli/Downloads/Theron/Include/ -I/opt/local/include/ boxspline2_parallel.cpp -L /opt/local/lib -lboost_thread-mt /Users/peterli/Downloads/Theron/Lib/libtheron.a
//
// Linux:
// mex -v -I/snle/lab/arch/Linux/lib/Theron/Include/ -I/usr/include/ boxspline2_parallel.cpp /Users/peterli/Downloads/Theron/Lib/libtheron.a /opt/matlab/bin/glnxa64/libboost_thread.so


#include "mex.h"
#include "math.h"
        
#define THERON_USE_BOOST_THREADS 1
#include "Theron/Framework.h"
#include "Theron/Actor.h"
#include "Theron/Receiver.h"
        
        
int DEFAULT_NUM_THREADS = 4;


struct SplineRequestMessage
{
  int start;
  int end;
  double *x;
  double *y;
  double *out;
};


class BoxSpliningActor : public Theron::Actor 
{    
public:
  inline BoxSpliningActor() {
    RegisterHandler(this, &BoxSpliningActor::SplineRequestHandler);
  }
  
private:
  inline void SplineRequestHandler(const SplineRequestMessage &message, const Theron::Address from) {
    for (int i = message.start; i < message.end; i++) {
      message.out[i] = boxspline2(message.x[i], message.y[i]);
    }

    TailSend(0, from);
  }
  
  inline float boxspline2(float x, float y) {
    float u = fabs(x) - fabs(y)/sqrt(3.0);
    float v = fabs(x) + fabs(y)/sqrt(3.0);
    
    if (u < 0) {
      u = -u;
      v = v + u;
    }
    if (2*u < v) u = v - u;
    
    if (v > 2.0) return 0.0; // Outside the support
    
    float g = u - v/2.0;
    if (v < 1.0) return 0.5 + ((5/3.0 - v/8.0)*v - 3)*v*v/4.0 + ((1 - v/4.0)*v + g*g/6.0 - 1)*g*g;
    if (u > 1.0) return (v - 2) * (v - 2) * (v - 2) * (g - 1) / 6.0;
    return 5/6.0 + ((1 + (1/3.0 - v/8.0)*v)*v/4.0 - 1)*v + ((1 - v/4.0)*v + g*g/6.0 - 1)*g*g;
  }
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check inputs
  if (nrhs < 2 || nrhs > 3)
    mexErrMsgIdAndTxt("boxspline2:nrhs", "Arguments should be the x and y coordinate matrices, plus optional number of threads");
  if (!mxIsDouble(prhs[0]))
    mexErrMsgIdAndTxt("boxspline2:prhs", "First argument must be a double matrix.");
  if (!mxIsDouble(prhs[1]))
    mexErrMsgIdAndTxt("boxspline2:prhs", "Second argument must be a double matrix.");
  
  // Check that x and y have same number of elements
  const mwSize numelx = mxGetNumberOfElements(prhs[0]);
  const mwSize numely = mxGetNumberOfElements(prhs[1]);
  if (numelx != numely)
    mexErrMsgIdAndTxt("boxspline2:prhs", "First and second arguments must have same number of arguments");

  int nthreads = DEFAULT_NUM_THREADS;
  if (nrhs > 2) {
    if (!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) > 1)
      mexErrMsgIdAndTxt("boxspline2:prhs", "Third argument should be scalar integer number of threads.");
    nthreads = mxGetScalar(prhs[2]);
  }

  double *inx = mxGetPr(prhs[0]);
  double *iny = mxGetPr(prhs[1]);
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
  double *out = mxGetPr(plhs[0]);

  Theron::Framework theron(nthreads);
  Theron::Receiver receiver;

  Theron::ActorRef spliners[nthreads];
  int step = numelx / nthreads;
  for (int i = 0; i < nthreads; i++) {
    SplineRequestMessage message;
    message.x = inx;
    message.y = iny;
    message.out = out;

    message.start = step * i;
    message.end = message.start + step;
    if (i == nthreads-1) message.end = numelx;

    spliners[i] = Theron::ActorRef(theron.CreateActor<BoxSpliningActor>());
    spliners[i].Push(message, receiver.GetAddress());
  }

  for (int i = 0; i < nthreads; i++) receiver.Wait();
}
