
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

function [f g H lcifs] = ll_func2_multi(p,basepars,stimpars,trainpars,kx_opt)

1;
nsets = size(stimpars.x,3);
f = 0;
g = zeros(size(p));
if (~(strcmp(basepars.filtermode,'nonsep') && strcmp(basepars.hessmode,'mult')))
    H = zeros(length(p));
else
    H = struct('kernel',zeros(basepars.fac*basepars.maxt,nsets),'kernel2',zeros(basepars.fac*basepars.maxt,nsets));
end
1;
for j=1:nsets % Add upt he value,gradient and hessians of the log-likelihood fns for each data set
    
   stimparsj = struct('x',stimpars.x(:,:,j),'dt',stimpars.dt);
   trainparsj.baseneuron_idx = trainpars.baseneuron_idx;
   trainparsj.analog_spikes = trainpars.analog_spikes;
   trainparsj.dt = trainpars.dt;
   trainparsj.psbasisGrad = trainpars.psbasisGrad{j};
   if (isfield(trainpars,'cpbasisGrad') && ~isempty(trainpars.cpbasisGrad))
       trainparsj.cpbasisGrad = trainpars.cpbasisGrad{j};
   else
       trainparsj.cpbasisGrad = [];
   end
   if (isfield(trainpars,'lgrad'))
       trainparsj.lgrad = trainpars.lgrad{j};
   end
   [sp_times trainparsj.D trainparsj.negSpikes] = pull_spikes(trainpars.analog_spikes,[],1,basepars.fac*basepars.maxt,stimpars.dt*basepars.maxt,trainpars.dt,(j-1)*stimpars.trialtime);        
   
   
   basepars.analog_frame_offset = (j-1)*stimpars.trialtime;
   
   
   1;
   
   [fj gj Hj lcifsj] = ll_func2(p,basepars,stimparsj,trainparsj);
   f = f + fj; g = g + gj;
   
   if (isnumeric(H))
      H = H + Hj;
   else
       H.kernel(:,j) = Hj.kernel;
       H.kernel2(:,j) = Hj.kernel2;
   end
   
end