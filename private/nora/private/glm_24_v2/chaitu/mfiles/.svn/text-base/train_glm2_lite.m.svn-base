% Training the parameters of the glm to data given a fixed set of stimulus filters

% Parameters to train:
% 1. stimulus filter (space x time) for each neuron
% 2. the coupling filters h (for each pair of neurons) - this is represented as cooridinates wrt some finite basis (of exponentials).
%    This is a 3D matrix (basis-coordinate-vector x neuron 1 x neuron 2)
% 3. the base firing rate of each neuron. This is a vector (neurons x 1)

% The function to be optimized is the log-likelihood of the data:

% sum_i[ sum_{t_{i,j}}[ log(L_{i,t})] - sum_t[L_{i,t}dt] ]


% the arguments:
% p0 - initial search point - the vector is blocked according to neurons.
%    - within each neuron block it is ordered as:
%    - base firing rate, stimulus filter, postspike basis coeffs, coupling basis coeffs (ordered by neuron)

% basepars - structure with the following fields:
    % maxt - # of timebins
    % n - spatial dimension
    % Nneurons - number of neurons 
    % postspike_basis - matrix where each column is a postspike basis function
    % couping_basis - "" coupling basis fns
    % nofilters_postspike - number of postspike basis fns
    % nofilers_coupling - "" coupling
    % fac - ratio of stimulus temporal resolution to spike train resolution
    %       - IMPORTANT!
    % fratio - ratio of timespan of postspike filter to stimulus filter
% stimpars - structure with the following fields:
    % x - stimulus (space x time)
    % dt - temporal resolution of stimulus
% trainpars - structure with the following fields:
    % dt - temporal resolution of spike train time bins
    % D - sparse logical matrix where ith column is the spike train of the ith neuron
% gopts - optimization options to pass to fminunc

% Returns:
% pstar - optimal solution returned by fminunc
% fstar, gstar,Hstar - function, gradient, and hessian value at the optimal point
% eflag,output - exit flag and output structure returned by fminunc

function [pstar fstar gstar Hstar eflag output] = train_glm2_lite(p0,basepars,stimpars,trainpars,gopts,NITER)

%-- Initialize the trainpars structure
% It contains any values that need not be computed on every optimization
% iteration such as the hessian spacetime constant matrix (i.e., the spike trains, temporally
% filtered stimulus), the convolutions of the spike trains with the
% postspike/coupling basis functions

if (~isfield(trainpars,'psbasisGrad'))
   fprintf('Computing ps basis grad...');
   trainpars.psbasisGrad = grad_basis(trainpars.D(:,1:basepars.Nneurons),basepars.postspike_basis);
   fprintf('done.\n');
end

if (size(trainpars.D,2) > 1 && ~isfield(trainpars,'cpbasisGrad'))
   fprintf('Computing cp basis grad...');
   trainpars.cpbasisGrad = grad_basis(trainpars.D,basepars.coupling_basis);
   fprintf('done.\n');
elseif (size(trainpars.D) == 1)
   trainpars.cpbasisGrad = [];
end

if (strcmp(basepars.filtermode,'nonsep') && ~isfield(trainpars,'lgrad'))
   fprintf('Computing lgrad for nonsep mode...');
   trainpars.lgrad = nonsep_lgrad(basepars,stimpars,trainpars); % edoi
   fprintf('done.\n');
end

%-- set up display during iteration
if (isempty(optimget(gopts,'PlotFcns'))) % plotting settings if not already set
    pfcn = setup_plotfcns(basepars.filtermode,basepars,stimpars,trainpars);
    gopts = optimset(gopts,'PlotFcns',pfcn);
end

%=== main part of training ===%
if (strcmp(basepars.hessmode,'mult') && strcmp(basepars.filtermode,'nonsep') && (~isfield(basepars,'normconstr') || basepars.normconstr == 0))
   % this is the case of my "nonsep" setting (edoi)
   gopts = optimset(gopts,'HessMult',@(Hinfo,Y) ll_hess_mult(Hinfo,Y,basepars,stimpars,trainpars));
   pstar = p0;
   for j=1:NITER
      [pstar, ~, eflag, output] = fminunc(@(p)ll_func2(p,basepars,stimpars,trainpars),pstar,gopts);
   end
   [fstar gstar] = ll_func2_lite(pstar,basepars,stimpars,trainpars);
   Hstar = [];%ll_hess_mult(Hinfo,eye(size(gstar,1)),basepars,stimpars,trainpars);
else
   fprintf('now running with lite version.  need to use train_glm2.m instead (edoi).\n')
end
