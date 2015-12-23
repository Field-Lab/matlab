% It returns function, gradient (and Hessian) value of the log likelihood
% of the data. This function can be used by fmincon

% pars has the following fields:
%
% Basic stuff:
% t - time vector
% maxt, n, Nneurons - time, spatial, neuronal dimensions
% x - stimulus
% D - data (3D binary matrix: time x trial x neuron)
% M - memory parameter

% History filter stuff:
% nofilters_postspike/coupling - number of postspike and coupling filter basis fns, resp.
% phi.postspike/coupling - phi params
% psi
% postspike/coupling_basis - matrix of postspike basis fns and coupling basis fns.

function [f g H lcifs] = ll_func2(p,basepars,stimpars,trainpars,kx_opt)

fl_filt_mode = ( strcmp(basepars.filtermode,'nonsep') || strcmp(basepars.filtermode,'raw'));

hessFlag = ~strcmp(basepars.hessmode,'mult') || ~fl_filt_mode;
% ~(1) || ~(1) = 0

% If any parameters are frozen at designated values, reset them to their designated values
frozenflag = isfield(basepars,'frozen_idx') && isfield(basepars,'frozen_vals');
if (frozenflag)
   p(basepars.frozen_idx) = basepars.frozen_vals;
end

% Calculate the log likelihood of the data
if (exist('kx_opt','var') && ~isempty(kx_opt))
   %fprintf('ll_func2: Supplying prefiltered stimulus...\n');
   [logprob lcifs cifs] = train_ll3(p,basepars,stimpars,trainpars,kx_opt);
else
   [logprob lcifs cifs] = train_ll3(p,basepars,stimpars,trainpars);
end

f = full(-logprob); % remember to negate - we want to maximize log likelihood
% Calculate gradient
if (nargout < 3)
   [g] = train_ll_grad4(p,basepars,stimpars,trainpars,lcifs,cifs);
else
   [g H] = train_ll_grad4(p,basepars,stimpars,trainpars,lcifs,cifs);
   if (hessFlag)
      H = -H;
   end
   %fprintf('*****')
   %keyboard, edoi
   %fprintf('Condition number of total hessian is %f\n',cond(H));
end

% If any parameters are frozen, set the gradient/hessian terms wrt these parameters to be zero.
if (frozenflag)
   g(basepars.frozen_idx) = 0;
   if (nargout >= 3 && hessFlag)
      H(basepars.frozen_idx,:) = 0;
      H(:,basepars.frozen_idx) = 0;
   end
end
g = -g;

% Determine if any stimulus filter regularization is specified
regflag = isfield(basepars,'stimreg') && strcmp(basepars.stimreg,'L2');
if (regflag)
   fprintf('## regularization is included ##\n')
   switch(basepars.filtermode)
      case {'sep_raw','rk2'}
         [f g H] = regularize_sep_spatial(p,basepars,f,g,H);    % Spatial filter only
         %[f g H] = regularize_sep_whole(p,basepars,f1,g1,H1)    % Whole sep filter
      case {'nonsep','raw'}
         [f g] = regularize_nonsep(p,basepars,f,g);
   end
end

%fprintf('f=%f    gnorm=%f\n',f,norm(g));
%[s1 t1 s2 t2] = get_sep_filters(p,basepars.n,basepars.Mk,[],[],0);
%fprintf('s1 norm = %0.5f , s2 norm = %0.5f , t1 norm = %0.5f , t2 norm = %0.5f\n',norm(s1),norm(s2),norm(t1),norm(t2));
