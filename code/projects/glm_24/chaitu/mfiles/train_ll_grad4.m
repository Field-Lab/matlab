% Function that computes the GRADIENT of the log likelihood of the data D for the
% parameters p

% pars has the following fields:
%
% Basic stuff:
% t - time vector
% maxt, n, Nneurons - time, spatial, neuronal dimensions
% kx - filtered stimulus (time x neuron)
% D - data (3D binary matrix: time x trial x neuron)
% M - memory parameter

% History filter stuff:
% nofilters_postspike/coupling - number of postspike and coupling filter basis fns, resp.
% phi.postspike/coupling - phi params
% psi
% postspike/coupling_basis - matrix of postspike basis fns and coupling basis fns.

% cifs are the corresponding rate fns for this data and given parameters -
% they will be used to calculate gradients

function [ll_grad ll_hess] = train_ll_grad4(p,basepars,stimpars,trainpars,lcifs,cifs)

% Basic variables
N = basepars.Nneurons;
Neff = size(trainpars.D,2);
ll_grad = zeros(size(p)); % this is the gradient vector

% Compute 'kernels' which are used as multipliers for computing
% gradients/hessian values. This is the factor that depends on the trial
% responses/cifs. It is multiplied by gradients/hessians of the input term
% (before the nonlinearity is applied).

fl_filt_mode = ( strcmp(basepars.filtermode,'nonsep') || strcmp(basepars.filtermode,'raw'));
% 1 (edoi)

hessFlag = (~fl_filt_mode || ~strcmp(basepars.hessmode,'mult') );
% = (0 || 0) = 0;

fprime = basepars.Nprime(lcifs);
kernel = fprime./cifs.*full(double(trainpars.D(:,trainpars.baseneuron_idx))) - trainpars.dt.*fprime; % approx. for dt small
if (nargout > 1)
    fdprime = basepars.Ndoubleprime(lcifs);
    kernel2 = (fdprime.*cifs - fprime.^2)./(cifs.^2).*full(double(trainpars.D(:,trainpars.baseneuron_idx))) - fdprime.*trainpars.dt;
    if (hessFlag)
        ll_hess = zeros(length(p),length(p));
    end
else
    kernel2 = [];
end

numpars = get_npars(basepars,Neff);
T = size(trainpars.D,1); % time dim of the spike train (fine resolution!)

for n=1:N % Compute gradient/hessian of parameters corresponding to each neuron (all cross terms in Hessian are 0)
    
    idx = ((n-1)*numpars+1):(n*numpars); % indices for this neuron in the param vector p
    
    % Make a numpars x pars.maxt matrix where the jth column is the
    % gradient of the input term at time j for this trial (for the nth neuron)
    
    % Compute the gradients of the input terms wrt the filters: (2*(nspace+pars.Mk) x pars.maxt)
    n_idx = trainpars.baseneuron_idx(n); % index of this neuron's spike train in trainpars.D
    neighbor_idx = [1:n_idx-1 n_idx+1:Neff]; % indices of the neighbor neuron's spike trains in trainpars.D
    
    if (Neff > 1)
        
        cpGrad = zeros(basepars.nofilters_coupling*(Neff-1),T);
        counter = 1;
        for n2 = neighbor_idx % neuron n is the one being AFFECTED
            cpGrad((counter-1)*basepars.nofilters_coupling+1:counter*basepars.nofilters_coupling,:) = trainpars.cpbasisGrad{n2}; % Add this neuron n2's spike train convolution
            counter = counter + 1;
        end
        
    else
        cpGrad = [];
    end
    
    
    if (isfield(basepars,'frozen_idx'))
        frozen_idx = basepars.frozen_idx;
    else
        frozen_idx = [];
    end

    1;
    
    % Compute the linear gradients if not in nonsep mode
    %if (~strcmp(basepars.filtermode,'nonsep'))
    if (~fl_filt_mode)
       %fprintf('**train_ll_grad4**') edoi
       %keyboard
        p_offset = get_pars_idx(basepars,n,Neff,'k');
        trainpars.lgrad{1} = linearfilt_grad(p,p_offset(1)-1,basepars,stimpars,trainpars,n);
        1;
        if (isfield(basepars,'XsqK') && basepars.XsqK)
            p_offset = get_pars_idx(basepars,n,Neff,'ksq');
            bp2 = basepars; bp2.padval = 0;%basepars.padval^2;
            %trainpars.lsqgrad{1} = linearfilt_grad(p,p_offset(1)-1,bp2,struct('dt',stimpars.dt,'x',(stimpars.x).^2),trainpars,n);
            trainpars.lsqgrad{1} = linearfilt_grad(p,p_offset(1)-1,bp2,struct('dt',stimpars.dt,'x',(sqrt((stimpars.x-0.5).^2))),trainpars,n);
        end
    elseif ((fl_filt_mode) && (isfield(basepars,'XsqK') && basepars.XsqK) )
            %trainpars.lsqgrad{1} = stimpars.dt .*(1/stimpars.dt.*trainpars.lgrad{1}).^2;
            trainpars.lsqgrad{1} = stimpars.dt .* sqrt(((1/stimpars.dt.*trainpars.lgrad{1})-0.5).^2);            
    end
    
    % Compute any extrinsic signal gradients, if necessary
    if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
       trainpars.extgrad{1} = interpmtx(basepars.interpfilt,basepars.ext_timepts,basepars.maxt)'; %sincinterpmtx(makeaxis(stimpars.dt,basepars.maxt),basepars.ext_timepts)';
    end
    
    
    1;
    if (nargout == 1 || (nargout>1 && ~hessFlag))
        [ll_grad(idx)] = multisample_gradHess_structs(basepars,n,kernel,[],trainpars,cpGrad,frozen_idx);
    elseif (hessFlag)
        1;
        [ll_grad(idx) ll_hess(idx,idx)] = multisample_gradHess_structs(basepars,n,kernel,kernel2,trainpars,cpGrad);
        if (~fl_filt_mode) % For separable filter case only, there are additional crossterms in the hessian
            ctoffset = get_pars_idx(basepars,n,Neff,'k');
            ll_hess = add_hess_crossterms(ll_hess,ctoffset(1)-1,basepars,stimpars,kernel,n);
            if (isfield(basepars,'XsqK') && basepars.XsqK)
                [nspace ntime] = get_nspacetime(basepars);
                ctoffset = get_pars_idx(basepars,n,Neff,'ksq');
                %ll_hess = add_hess_crossterms(ll_hess,ctoffset(1)-1,basepars,struct('dt',stimpars.dt,'x',(stimpars.x).^2),kernel,n);
                ll_hess = add_hess_crossterms(ll_hess,ctoffset(1)-1,basepars,struct('dt',stimpars.dt,'x',sqrt((stimpars.x-0.5).^2)),kernel,n);                
            end
        end
    end
end

if (nargout >1 && ~hessFlag)
    % fprintf('Creating Hinfo structure\n');
    % note: because ll_hess is structure instead of matrix, you'll get some
    % error if 'funvalcheck' is on.  (edoi)
    ll_hess = struct('kernel',kernel,'kernel2',kernel2);
end

