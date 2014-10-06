#include <iostream>
#include "mex.h"
#include "cppa/cppa.hpp"
using namespace cppa;


void slogger() {
  receive(on<int>() >> [](int value) {
    reply((value * 20) + 2);
  });
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto sl = spawn(slogger);
  send(sl, 2);
  receive(on<int>() >> [](int value) {
    mexPrintf("%d\n", value);
  });
}