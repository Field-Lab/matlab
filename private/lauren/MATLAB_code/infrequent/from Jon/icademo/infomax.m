function [signals, A, W] = infomax(mixes)
% INFOMAX     Performs ICA using Bell-Sejnowski "Info-Max" algorithm
%
% Usage:   [signals, A, W] = infomax(mixes)
%                [signals] = infomax(mixes)
%
% Inputs:     mixes  -  NxP matrix with N sources and P data
%                          points in each source.
%
% Outputs:        A  -  mixing matrix (recovered)    [NxN matrix]
%                 W  -  unmixing matrix (recovered)  [NxN matrix]
%           signals  -  unmixed sources              [NxP matrix]
%
% See also:  NATGRAD
%
% INFOMAX performs Indepedent Component Analysis (ICA) on a data
% set. The goal of ICA is to find the linear transformation which
% projects the data to a new coordinate system that maximizes the
% statistical independence between axes. These axes do not (and are
% often not) orthogonal! A great introduction is presented in:
%     Bell & Sejnowski (1995) Neural Computation 7(6) 1129-1159.
%
% This algorithm presented in the article above was subsequently
% refined, sped up and made more biologically plausible (for
% applications to neuroscience) in:
%     Amari et al (1996) Advances in NIPS (Volume 8)
%
% The current algorithm employed in this function is fully detailed
% in the above two papers. A few refinements and a clear derivation
% is presented in:
%     Bell & Sejnowski (1997) Vision Research 37(23) 3327-3338.
% The final algorithm presented here uses the equations in the
% above paper. In addition, the variable names are consistent with
% the paper as well.
%
%    W  - weight matrix                    X  - input matrix
%    Wz - sphering matrix                  U  - outputs (U = W * X)
%    A  - mixing matrix ( A=inv(W * Wz) )  Y  - sigmoid(U)  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A Version History
% -----------------
%
% Original code by Tony Bell is available at:
%      ftp://ftp.cnl.salk.edu/pub/tony/sep96.public
%
% Code adapted and (slightly) modified.
%      23 August 2002,  Jon Shlens  |  jonshlens [at] ucsd.edu
%


%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. FIND SIZE OF SOURCES
%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,P] = size(mixes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. PERMUTE SOURCES OVER TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If the signals are time series, then one needs to assure
%  stationarity in the data set. Shuffling the time points of all
%  signals provides a time-scrambled version of the inputs.

% generate the permutation matrix (across columns)
permute = randperm(P);

% permute the data matrix
X = mixes(:,permute);


%%%%%%%%%%%%%%%%%%%%
% 3. SPHERE THE DATA
%%%%%%%%%%%%%%%%%%%%
%  Pre-whitening the data ensures that the learned unmixing matrix
%  W is an orthonormal, rotation matrix with no scaling factor.
%  Mathematically, cov(X')=4*eye(N). The reason for pre-processing
%  the data in this fashion is because empirically the sigmoid
%  nonlinearity usually settles on creating a U matrix (where U=WX)
%  where cov(U')=4I. For convergence speed we want the infomax
%  algorithm to solely focus on finding the proper rotation without
%  the scaling factors (due to variance). Hence, the "sphereing"
%  step manually performs the scaling change which the sigmoid 
%  nonlinearity would have settled on anyways.

% subtract off the mean of each mixed source
X = X - repmat(mean(X,2),1,P);

% calculate the decorrelating matrix
% Note: When trying the math out, remember the covariance matrix
% is symmetric.
Wz = 2 * inv(sqrtm(cov(X')));

% whiten X so that: cov(X') = 4*eye(N)
X = Wz * X;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. INITIALIZE THE MIXING MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  For this form of ICA, the number of outputs equals the number of
%  inputs. Therefore, this unmixing matrix will be square. We can
%  choose a random matrix (i.e. rand(N)), however, the identity
%  matrix will suffice since the data is pre-whitened.

W = eye(N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. RUN INFOMAX LEARNING ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Use annealing method that is common to Connectionist
%  learning. Basically, start with a fast learning rate and settle
%  down over successive runs.

% -- Learning Rate: the i'th value for the i'th sweep
learning_rate = [ 100 10 ones(1,5000) ] * 1e-4;

% Perform multiple learning sweeps of the data until it hopefully
% converges. In practice, this can be more of an art than a
% science.

for l = learning_rate

  W = natgrad(W, X, l);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. CALCULATE OUTPUT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Un-mixing matrix  (recovered)
W = W * Wz;

% The mixing matrix (recovered)
A = inv(W);

% Un-mixed signals
signals = W * mixes;
