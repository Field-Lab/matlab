% It returns function, gradient (and Hessian?) value of the log likelihood
% of the data. This function can be used by fmincon

% pars has the following fields:
%
% Basic stuff:
% t - time vector
% maxt,n,Nneurons - time, spatial, neuronal dimensions
% x - stimulus
% D - data (3D binary matrix: time x trial x neuron)
% M - memory parameter

% History filter stuff: 
% nofilters_postspike/coupling - number of postspike and coupling filter basis fns, resp.
% phi.postspike/coupling - phi params
% psi
% postspike/coupling_basis - matrix of postspike basis fns and coupling basis fns.

function [f g cifs] = ll_func_cg(p,basepars,stimpars,trainpars,kx_opt)

% Calculate the log likelihood of the data
if (exist('kx_opt','var') && ~isempty(kx_opt))
    %fprintf('ll_func2: Supplying prefiltered stimulus...\n');
    [logprob cifs] = train_ll3(p,basepars,stimpars,trainpars,kx_opt);
else
    [logprob cifs] = train_ll3(p,basepars,stimpars,trainpars);
end

f = full(-logprob); % remember to negate - we want to maximize log likelihood
[g] = train_ll_grad_hessmult(p,basepars,stimpars,trainpars,cifs);
g = -g;
%fprintf('f=%f    gnorm=%f\n',f,norm(g));





