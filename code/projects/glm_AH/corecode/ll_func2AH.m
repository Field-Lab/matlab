% It returns function, gradient (and Hessian) value of the log likelihood
% of the data. This function can be used by fmincon
%% CALLS TRAIN_LL3AH   AND TRAIN_LL_GRAD4


%%% sum_i[ sum_{t_{i,j}}[ log(L_{i,t})] - sum_t[L_{i,t}dt] ]   %% what is i
%%% in this case???

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

function [f g H lcifs] = ll_func2AH(p,Basepars,Stimpars,Trainpars,kx_opt)
%[f g H lcifs] = ll_func2AH(p,Basepars,Stimpars,Trainpars)


% 3 OR 4 SECONDS SPENT HERE    UNCOUPLED
% Calculate the log likelihood of the data given kx_opt exists
if (exist('kx_opt','var') && ~isempty(kx_opt))
   %fprintf('ll_func2: Supplying prefiltered stimulus...\n');
   [logprob lcifs cifs] = train_ll5_cpAH(p,Basepars,Stimpars,Trainpars,kx_opt);
else
   [logprob lcifs cifs] = train_ll5_cpAH(p,Basepars,Stimpars,Trainpars);  %%%%
end

% FULL IS MATLAB FUNCTION WHICH CHANGES SPARSE MATRIX TO A FULL MATRIX 
f = full(-logprob); % remember to negate - we want to maximize log likelihood

filtermode = ( strcmp(Basepars.k_filtermode,'nonsep') || strcmp(Basepars.k_filtermode,'raw'));
hessFlag = ~strcmp(Basepars.hessmode,'mult') || ~filtermode;  
%   THIS TAKES ONLY 15 SECONDS    UNCOUPLED 
if (nargout < 3)
   [g] = train_ll_grad5AH(p,Basepars,Stimpars,Trainpars,lcifs,cifs);
else
   [g H] = train_ll_grad5AH(p,Basepars,Stimpars,Trainpars,lcifs,cifs);
   if (hessFlag)
      H = -H;
   end
   
end




%% hessFlag is 1 for rank-2
%% something about flipping the hessian as well
% If any parameters are frozen at designated values, reset them to their designated values
frozenflag = isfield(Basepars,'frozen_idx') && isfield(Basepars,'frozen_vals');
if (frozenflag)
   p(Basepars.frozen_idx) = Basepars.frozen_vals;
end



% If any parameters are frozen, set the gradient/hessian terms wrt these parameters to be zero.
if (frozenflag)
   g(Basepars.frozen_idx) = 0;
   if (nargout >= 3 && hessFlag)
      H(Basepars.frozen_idx,:) = 0;
      H(:,Basepars.frozen_idx) = 0;
   end
end
g = -g;


% Determine if any stimulus filter regularization is specified
regflag = isfield(Basepars,'stimreg') && strcmp(Basepars.stimreg,'L2');
if (regflag)
   fprintf('## regularization is included ##\n')
   switch(Basepars.k_filtermode)
      case {'sep_raw','rk2'}
         [f g H] = regularize_sep_spatial(p,Basepars,f,g,H);    % Spatial filter only
         %[f g H] = regularize_sep_whole(p,Basepars,f1,g1,H1)    % Whole sep filter
      case {'nonsep','raw'}
         [f g] = regularize_nonsep(p,Basepars,f,g);
   end
end

%load(sprintf('%s/Track_Progress.mat',Basepars.trackprog_dir))
%ind = Track_Progress(1).counter+1;
%Track_Progress(ind).p       = p;
%Track_Progress(ind).g       = g;
%Track_Progress(ind).H       = H;
%Track_Progress(ind).LogProb = logprob;
%Track_Progress(1).counter = ind;
%save(sprintf('%s/Track_Progress.mat',Basepars.trackprog_dir ),'Track_Progress');
%clear Track_Progress

end




%fprintf('f=%f    gnorm=%f\n',f,norm(g));
%[s1 t1 s2 t2] = get_sep_filters(p,Basepars.n,Basepars.Mk,[],[],0);
%fprintf('s1 norm = %0.5f , s2 norm = %0.5f , t1 norm = %0.5f , t2 norm = %0.5f\n',norm(s1),norm(s2),norm(t1),norm(t2));
