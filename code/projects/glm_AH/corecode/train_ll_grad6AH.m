%%% BASICALLY A WRAPPER FUNCTION %%%%%%

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

function [ll_grad ll_hess] = train_ll_grad6AH(p,Basepars,Stimpars,Trainpars,lcifs,cifs)



nNeighbors = length(Basepars.cp_Neighbors);
filtermode = ( strcmp(Basepars.k_filtermode,'nonsep') || strcmp(Basepars.k_filtermode,'raw') );
hessFlag = (~filtermode || ~strcmp(Basepars.hessmode,'mult') );




%%% KERNEL_GRAD AND KERNEL_HESS %%%%%%%%%%%
% (MICROBINS,1) GENERAL VECTOR THAT DOTS INTO PARAM SPECIFIC DERIVATIVES
ll_grad     = zeros(size(p));
kernel_grad = full(double(Trainpars.logicalspike_microbin_Home)) - Trainpars.dt*cifs; % approx. for dt small
if (nargout > 1)
    kernel_hess = -Trainpars.dt *cifs;
    if (hessFlag)
        ll_hess = zeros(length(p),length(p));
    end
else
    kernel_hess = [];
end
numpars = length(Basepars.p0);
microbins = size(Trainpars.logicalspike_microbin_Home,1); 
idx = 1:numpars;    
if (isfield(Basepars,'frozen_idx'))
	frozen_idx = Basepars.frozen_idx;
else
	frozen_idx = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% LINEARFILT_GRAD   %%%%%%%%%%%%
%%%% For rk2 [Xkt;Xks;Xkt;Xks] %%%%%%
if strcmp( Basepars.k_filtermode ,'rk2' )
    if (~filtermode)
        %%%%%%%  old version  stimpars.dt * [Xkt;Xks:Xkt;Xks]  
        %%%%%%%  new version                [Xkt;Xks;Xkt;Xks]
        Trainpars.lGrad = linearfilt_grad2AH(p,Basepars,Stimpars,Trainpars,1);
        %%%% TRY WRITING A NEW LINEARFILT_GRADAH2   WITHOUT THE EXTRA
        %%%% STIM.DT!!!! %%%%

        if (isfield(Basepars,'XsqK') && Basepars.XsqK)
            p_offset = get_pars_idx(Basepars,n,Neff,'ksq');
            bp2 = Basepars; bp2.padval = 0;%Basepars.padval^2;
                %Trainpars.lsqgrad{1} = linearfilt_grad(p,p_offset(1)-1,bp2,struct('dt',Stimpars.dt,'x',(Stimpars.movie_ROI).^2),Trainpars,n);
            Trainpars.lsqgrad = linearfilt_grad2AH(p,p_offset-1,bp2,struct('dt',Stimpars.dt,'x',(sqrt((Stimpars.movie_ROI-0.5).^2))),Trainpars,1);
        end
    elseif ((filtermode) && (isfield(Basepars,'XsqK') && Basepars.XsqK) )
        Trainpars.lsqgrad = Stimpars.dt .* sqrt(((1/Stimpars.dt.*Trainpars.lgrad{1})-0.5).^2); 
    end
elseif strcmp( Basepars.k_filtermode ,'STA' )
    Trainpars.lGrad= Basepars.kx_STA';
end

if isfield(Basepars, 'rect_conv_spSTA')
    Trainpars.recGrad = Basepars.rect_conv_spSTA';
end

if isfield(Basepars, 'rect_full')
    Trainpars.recGrad = Basepars.rect_full;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  
% Compute any extrinsic signal gradients, if necessary %
if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts))
       Trainpars.extgrad{1} = interpmtx(Basepars.interpfilt,Basepars.ext_timepts,Basepars.maxt)'; %sincinterpmtx(makeaxis(Stimpars.dt,Basepars.maxt),Basepars.ext_timepts)';
end

if (nargout == 1 || (nargout>1 && ~hessFlag))
    [ll_grad(idx)] = multisample_gradHess_structsAH_rec(Basepars,Stimpars,Trainpars,kernel_grad,kernel_hess,frozen_idx);
elseif (hessFlag) 
	[ll_grad(idx) ll_hess(idx,idx)] = multisample_gradHess_structsAH_rec(Basepars,Stimpars,Trainpars,kernel_grad,kernel_hess,frozen_idx);
end

if (nargout >1 && ~hessFlag)
    % fprintf('Creating Hinfo structure\n');
    % note: because ll_hess is structure instead of matrix, you'll get some
    % error if 'funvalcheck' is on.  (edoi)
    ll_hess = struct('kernel_grad',kernel_grad,'kernel_hess',kernel_hess);
end




%%% UNNECESARRY AND COSTLY USE OF COMPUTING TIME 
%%% DO NOT USE THIS OPERTATION!! 
%{

if (nNeighbors > 0)      
	cpGrad = zeros(Basepars.cp_filternumber*(nNeighbors),microbins);
    counter = 1;
	for iNeighbor = 1:nNeighbors % neuron n is the one being AFFECTED
        tic
        cpGrad((counter-1)*Basepars.cp_filternumber+1:counter*Basepars.cp_filternumber,:) = Trainpars.cpbasisGrad{iNeighbor}; % Add this neuron n2's spike train convolution
        counter = counter + 1;
        toc
	end
else
	cpGrad = [];
end
%}    
%if (~filtermode) % For separable filter case only, there are additional crossterms in the hessian
%    ctoffset = get_pars_idx2AH(Basepars,'mu');
%    ll_hess = add_hess_crosstermsAH(ll_hess,ctoffset(1)-1,Basepars,Stimpars,kernel_grad,n);
%    if (isfield(Basepars,'XsqK') && Basepars.XsqK)
%        [nspace ntime] = get_nspacetimeAH(Basepars);
%       % ctoffset = get_pars_idxAH(Basepars,n,Neff,'ksq');
%        %ll_hess = add_hess_crossterms(ll_hess,ctoffset(1)-1,Basepars,struct('dt',Stimpars.dt,'x',(Stimpars.movie_ROI).^2),kernel,n);
%        ll_hess = add_hess_crosstermsAH(ll_hess,ctoffset(1)-1,Basepars,struct('dt',Stimpars.dt,'x',sqrt((Stimpars.movie_ROI-0.5).^2)),kernel_grad,n);                
%    end
%end

end

