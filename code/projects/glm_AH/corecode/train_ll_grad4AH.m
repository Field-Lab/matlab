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

function [ll_grad ll_hess] = train_ll_grad4AH(p,Basepars,Stimpars,Trainpars,lcifs,cifs)

nNeighbors = length(Basepars.cp_Neighbors);
ll_grad = zeros(size(p)); % this is the gradient vector

% Compute 'kernels' which are used as multipliers for computing
% gradients/hessian values. This is the factor that depends on the trial
% responses/cifs. It is multiplied by gradients/hessians of the input term
% (before the nonlinearity is applied).

% DETERMINES IF WE NEED EXTRA CROSS TERMS FOR K_FILTER
filtermode = ( strcmp(Basepars.k_filtermode,'nonsep') || strcmp(Basepars.k_filtermode,'raw'));

hessFlag = (~filtermode || ~strcmp(Basepars.hessmode,'mult') );
% = (0 || 0) = 0;

fprime = Basepars.Nprime(lcifs);
kernel = fprime./cifs.*full(double(Trainpars.logicalspike_microbin_Home(:,Trainpars.baseneuron_idx))) - Trainpars.dt.*fprime; % approx. for dt small
if (nargout > 1)
    fdprime = Basepars.Ndoubleprime(lcifs);
    kernel2 = (fdprime.*cifs - fprime.^2)./(cifs.^2).*full(double(Trainpars.logicalspike_microbin_Home(:,Trainpars.baseneuron_idx))) - fdprime.*Trainpars.dt;
    if (hessFlag)
        ll_hess = zeros(length(p),length(p));
    end
else
    kernel2 = [];
end

numpars = get_nparsAH(Basepars,Neff);
T = size(Trainpars.logicalspike_microbin_Home,1); % time dim of the spike train (fine resolution!)

for n=1:N % Compute gradient/hessian of parameters corresponding to each neuron (all cross terms in Hessian are 0)
    
    idx = ((n-1)*numpars+1):(n*numpars); % indices for this neuron in the param vector p
    
    % Make a numpars x pars.maxt matrix where the jth column is the
    % gradient of the input term at time j for this trial (for the nth neuron)
    
    % Compute the gradients of the input terms wrt the filters: (2*(nspace+pars.Mk) x pars.maxt)
    n_idx = Trainpars.baseneuron_idx(n); % index of this neuron's spike train in Trainpars.D
    neighbor_idx = [1:n_idx-1 n_idx+1:Neff]; % indices of the neighbor neuron's spike trains in Trainpars.logicalspike_microbin_Home
    
    if (Neff > 1)
        
        cpGrad = zeros(Basepars.nofilters_coupling*(Neff-1),T);
        counter = 1;
        for n2 = neighbor_idx % neuron n is the one being AFFECTED
            cpGrad((counter-1)*Basepars.nofilters_coupling+1:counter*Basepars.nofilters_coupling,:) = Trainpars.cpbasisGrad{n2}; % Add this neuron n2's spike train convolution
            counter = counter + 1;
        end
        
    else
        cpGrad = [];
    end
    
    
    if (isfield(Basepars,'frozen_idx'))
        frozen_idx = Basepars.frozen_idx;
    else
        frozen_idx = [];
    end

    1;
    
    % Compute the linear gradients if not in nonsep mode
    %if (~strcmp(Basepars.filtermode,'nonsep'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%% LOOK HERE TO UPDATE THE CODE FOR NEW PARAM
    %%%%%%%%%%%%%%%%%%%%%%%%%%% NAMES  %%%%%%%%%%%%%%%%%%%%%
    if (~filtermode)
       %fprintf('**train_ll_grad4**') edoi
       %keyboard
        p_offset = get_pars_idxAH(Basepars,n,Neff,'k');
        Trainpars.lgrad{1} = linearfilt_gradAH(p,p_offset(1)-1,Basepars,Stimpars,Trainpars,n);
        if (isfield(Basepars,'XsqK') && Basepars.XsqK)
            p_offset = get_pars_idx(Basepars,n,Neff,'ksq');
            bp2 = Basepars; bp2.padval = 0;%Basepars.padval^2;
            %Trainpars.lsqgrad{1} = linearfilt_grad(p,p_offset(1)-1,bp2,struct('dt',Stimpars.dt,'x',(Stimpars.movie_ROI).^2),Trainpars,n);
            Trainpars.lsqgrad{1} = linearfilt_grad(p,p_offset(1)-1,bp2,struct('dt',Stimpars.dt,'x',(sqrt((Stimpars.movie_ROI-0.5).^2))),Trainpars,n);
        end
    elseif ((filtermode) && (isfield(Basepars,'XsqK') && Basepars.XsqK) )
            %Trainpars.lsqgrad{1} = Stimpars.dt .*(1/Stimpars.dt.*Trainpars.lgrad{1}).^2;
            Trainpars.lsqgrad{1} = Stimpars.dt .* sqrt(((1/Stimpars.dt.*Trainpars.lgrad{1})-0.5).^2);            
    end
    
    % Compute any extrinsic signal gradients, if necessary
    if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts))
       Trainpars.extgrad{1} = interpmtx(Basepars.interpfilt,Basepars.ext_timepts,Basepars.maxt)'; %sincinterpmtx(makeaxis(Stimpars.dt,Basepars.maxt),Basepars.ext_timepts)';
    end
    
    
    1;
    if (nargout == 1 || (nargout>1 && ~hessFlag))
        % THIS IS IT IF IT IS NOT RANK 2
        [ll_grad(idx)] = multisample_gradHess_structsAH(Basepars,n,kernel,[],Trainpars,cpGrad,frozen_idx);
    elseif (hessFlag)
       
        % SYMMETRIC PORTION OF THE HESSIAN.. NO CROSS TERMS YET    ~15 secs
        [ll_grad(idx) ll_hess(idx,idx)] = multisample_gradHess_structsAH(Basepars,n,kernel,kernel2,Trainpars,cpGrad);
       
        if (~filtermode) % For separable filter case only, there are additional crossterms in the hessian
            ctoffset = get_pars_idxAH(Basepars,n,Neff,'k');
            ll_hess = add_hess_crosstermsAH(ll_hess,ctoffset(1)-1,Basepars,Stimpars,kernel,n);
            if (isfield(Basepars,'XsqK') && Basepars.XsqK)
                [nspace ntime] = get_nspacetimeAH(Basepars);
                ctoffset = get_pars_idxAH(Basepars,n,Neff,'ksq');
                %ll_hess = add_hess_crossterms(ll_hess,ctoffset(1)-1,Basepars,struct('dt',Stimpars.dt,'x',(Stimpars.movie_ROI).^2),kernel,n);
                ll_hess = add_hess_crosstermsAH(ll_hess,ctoffset(1)-1,Basepars,struct('dt',Stimpars.dt,'x',sqrt((Stimpars.movie_ROI-0.5).^2)),kernel,n);                
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

