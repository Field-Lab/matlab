function W = natgrad(W, X, lr)
% NATGRAD     Natural gradient learning rule for "Info-Max" ICA
%
% Usage:   W_new = natgrad(W_old, X, lr)
%
% Variables:    W_old, W_new  -  the new and old weight matrices
%                        X    -  the mixed sources
%                        lr   -  the learning rate
%
% See also:  INFOMAX
%
% Performs natural gradient ascent on the weight matrix W in order
% find the maximal joint entropy of the nonlinearity (sigmoid)
% applied to the W*X (i.e. sigmoid(W*X) ).
%
% NATGRAD uses batch learning by breaking up large data sets into
% "blocks." This fasciliates faster and more stable learning of the
% weight matrix W.
%
% Bell (the original author) suggested a learning rate of 0.01 for
% 2 source/signal separation. For higher dimensional data sets
% (i.e. >2), smaller learning rates are needed to avoid instability
% (although sacrificing convergence speed). Finally, regardless
% of the dimensionality, all learning should include some anneealing
% where by the learning rate is gradually decreased over time in order
% to fine tune the solution.
%
% Note: This learning rule makes use of the natural gradient trick
% (multiplying the gradient by W^T*W) discovered and detailed in:
%     Amari et al (1996) Advances in NIPS (Volume 8)
%
% This learning rule is detailed in:
%     Bell & Sejnowski (1997) Vision Research 37(23) 3327-3338.
%
% Original code by Tony Bell available at:
%      ftp://ftp.cnl.salk.edu/pub/tony/sep96.public
%
% Code adapted and slightly modified:
%      23 August 2002,  Jon Shlens  |  jonshlens [at] ucsd.edu
%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. SET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%

% general paramters
[N,P] = size(X);    % size of mixed signals
I = eye(N);         % identity matrix of [NxN]

% batch learning paramters (nothing fancy)
if P > 100
  %%%%%% -----> implement batch learning! 
  b = 30;           % block size (# of data pts) in batch
else
  %%%%%  -----> no batch learning.
  b = 1;            % b=1 reduces down to standard learning
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. A BIG NOTE ABOUT BATCH LEARNING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is just one commented note with no code. I feel it
% necessary because there is a big, subtle *trick* below in an
% equation which needs some explaining and is not detailed in any
% paper.
%
% Notice the calculation of "dW" has a term "b*I". If you are
% diligent reader, you might notice that in the 1997 paper this
% term should just be "I" (the identity matrix) and not multiplied
% by "b" (the number of data points in a batch.
%
% The not-so-obvious reason for this (which I had to grapple with
% for a while) is that this is *batch* learning. To see this, let
% me outline the *standard* learning from 1997 paper. In the standard 
% method, the algorithm first selects ONE (just 1) data point (all N 
% dimensions!) and calculates the associated dW. It then updates W.
% And, then continues on with another data point, etc.
%
% In batch learning the algorithm first selects ONE data point and
% calculates the associated dW. The algorithm then selects another
% data point (it does not update W!) and calculates dW. The
% algorithm does this for all "b" points (the number of points in a
% batch). After calculating all "b" "dW"'s, it then sums them up to
% form 1 dW. The weight matrix W is then updated with this summed
% dW. Good, eh?
%
% Now back to the weird equation. What does an algorithm for this
% batch learning look like?  Well, one could save all "b" "dW"'s,
% add them up, and then change W. This is exactly what the weird
% term (b*I) does!
%
%  Consider the case of b=1. This reduces it down to equation
%  13. This is the equation for 1 data point. Now let's say that
%  our batch size is 10 (b=10). This means we need to take 10
%  equations for dW and add them together. Well, let's do it. If we
%  write out the math, we see that the right-hand-side W pops out
%  (distributive property), and we are left with 10*I and a sum
%  over the Y_hat*U' for all 10 data points. The sum over Y_hat*U'
%  is what the second term is (matrix multiplication). And the
%  first term is 10*I. Thus, the "b" is the correct multiplicative
%  factor before "I" in batch learning.


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. RUN BATCH LEARNING
%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each "b" data points in the series
for t = 1 : b : fix(P/b)

  % find range of data points in X for batch
  first = t;
  last  = t + b-1;

  % calculate U in U = WX          (equation 2 in 1997 paper)
  U = W * X(:, first:last );

  % calculate Y_hat for sigmoid    (equation 12 in 1997 paper)
  Y_hat = 1 - 2./(1+exp(-U));

  % calculate batch weight change  (equation 13 in 1997 paper)
  % Note: see huge note above!
  dW = (b*I + Y_hat * U') * W;

  % apply weight change
  W = W + lr * dW;

end;
