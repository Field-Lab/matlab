% cleaned and altered again... 10-29-2012

%%% SAME AS TRAIN_LL4_CP EXCEPT NOW CALLS FILTERSTIMULUS_TRAIN2 
%%% INSTEAD OF FILTERSTIMULUS_TRAIN


% calls filterstimulus_train2AH


% Written on 10-14  AK Heitman   
% Tis version just uses spikes from train.cpbasisgrad   or
% train.psbasisgrad
%%%%%%%%% 
%%%%%%%%% CALLS ll_eval_mex %%%%%%%%%%%%
%%% but doesnt have to.. cause matlab code is just as fast %%


%INPUT_TERM MEANS lcifs FOR NOW

% FUNCTION CALLS: FILTERSTIMULUSTRAIN, GET_NPARS, GET_PARS_IDX
% ALSO CALLS GENINTERP (BUT NOT IN MY PARAM RANGE)


% Returns the log probability of the data, the CIFS, and the filtered
% stimulus (under current parameter settings).

% The stimulus and the spike train can have different resolutions
% the dt should be specified correctly in pars and Stimpars, resp.

% GENERATE  A KX A PS AND A CP  AND COMPUTE FROM THERE ?? 
% NEED TO KEEP TRACK OF INPUT_TERM .. BECOMES THE FINAL LCIFS
function [log_prob lcifs cifs kx] = train_ll5_cpAH(p,Basepars,Stimpars,Trainpars,kxopt)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(Basepars, 'LL_eval_C')
    Basepars.LL_eval_C = false ;
end
nNeighbors = length(Basepars.cp_Neighbors);
[spacePixels stimFrames] = size(Stimpars.movie_ROI);
numpars = length(Basepars.p0);
microbins = stimFrames * Basepars.spikebins_perstimframe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% KX + MU  TERM OF THE lcifs %%%%
mu_index = Basepars.paramind.MU;
mu = p(mu_index);% 1  x N   %% looks like the mu term of the filter
%%% Run through stimulus filter portion  (no post spike or coupling yet)
if (~exist('kxopt','var'))
    if (isfield(Basepars,'XsqK') && Basepars.XsqK)
        %fprintf('Using square filter\n');
        1;
        [kx kxsq] = filterstimulus_train2AH(p,Basepars,Stimpars); % T x N
        kx = kx + kxsq;
    else
        kx = filterstimulus_train2AH(p,Basepars,Stimpars);
    end
else
    Lind = Basepars.paramind.L;
   
    kx = p(Lind)*kxopt';
end
kxmu_stimblocked =( kx + repmat(mu,1,stimFrames) )'; % T x N   %%% AH i'm pretty sure this is analogous to the bas firing addition
kxmu_microbin = reprows(kxmu_stimblocked ,Basepars.spikebins_perstimframe); 
if(~isfield(Basepars,'frame_offset'))
    Basepars.frame_offset = 0;
end
if (isfield(Basepars,'analog_frame_offset') && ~isempty(Basepars.analog_frame_offset))
    numbinsadded = double(int32(Basepars.analog_frame_offset/Trainpars.dt));
else
    numbinsadded = Basepars.frame_offset*Basepars.spikebins_perstimframe;    
end
kxmu_microbin= [zeros(numbinsadded,size(kxmu_microbin,2)); kxmu_microbin];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% POST SPIKE FILTERING %%%%%%%%
if isfield(Basepars,'ps_FIXconvolved')
    psidx    = Basepars.paramind.PS; 
    PSupdate = p(psidx) * Basepars.ps_FIXconvolved;
else
    psidx = Basepars.paramind.PS;  %indexes the ps indices of p
    ps_weights = p(psidx(:));
    PS = Basepars.ps_basis*ps_weights; % Basepars.Mhist x N    
    PSupdate = (ps_weights') *cell2mat(Trainpars.psbasisGrad); 
end
%%% COUPLING HAPPENS HERE !!!  %%%%%
%keyboard
if (Basepars.Coupling && nNeighbors> 0)
    cpidx = Basepars.paramind.CP;
    cp_weights = reshape(p(cpidx(:)),Basepars.cp_filternumber,nNeighbors); % THIS NEEDS TO BE CONSISTENT IN WHAT P MAPS TO!! %%%
    CP = Basepars.cp_basis*cp_weights;  %% CP is a matrix !!
    
    CPupdate = zeros(microbins, 1);
    for iNeigh = 1 : nNeighbors        
        addterm  = cp_weights(:,iNeigh)' * (Trainpars.cpbasisGrad{iNeigh});
        CPupdate = CPupdate + addterm';
      %  figure; plot(addterm(1:500));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FINALIZE THE lcifs and CIF %%%%%%%%%%%
lcifs = kxmu_microbin + PSupdate';
if (Basepars.Coupling && nNeighbors> 0) 
    lcifs = lcifs + CPupdate;
end


%%% we should stay away from over a million htz type cifs..period
%{
toohigh = find ( lcifs > 100);
if ~isempty(toohigh)
    display('Change initial estimates .. CIF blows up.. look train_ll5 of llfunc2')
    lcifs(toohigh) = 100;
end

%%% don't want to get too close to 
toolow = find ( lcifs < -33 );
if ~isempty(toolow)
    display('Change initial estimates .. CIF approaches machine epsilon.. look train_ll5 ofl lfunc2')
    lcifs(toolow) = -25;
end
%}    
%%% by defn  duh! 
cifs = exp(lcifs);  % NEED TO BE CAREFUL OF THIS RUNNING OF TO INFINITY
%if (sum(cifs(:)<eps)>0)
%    cifs = max(eps,cifs);
%end
%if (sum(cifs(:)>(1/Trainpars.dt))>0)
%    cifs = min(1/Trainpars.dt,cifs);
%end


%%%% IDEA : CAP FIRING RATE AT 1000 HZ FOR
%%%% USE THIS TO CAP LCIFS
%%%% NEED BASEPARS.MAXFIRITNG RATE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  COMPUTING THE LOG PROB    %%%%%%%%%%%%%%
%%%%%%%%%%%%% THE FUNCTION WE OPTIMIZE   %%%%%%%%%%%%%%
spike_term =  sum( lcifs .* Trainpars.logicalspike_microbin_Home );
nospike_term = -Trainpars.dt* sum(cifs); 

if (abs(max(cifs(:))) > 10^6)
    fprintf('Warning: large max cif value: logC = %f. spike_term=%f nospike_term=%f\n',full(log(max(cifs(:)))),full(spike_term),full(nospike_term));
end
if (~isreal(nospike_term))
    fprintf('Warning: No spike is imaginary - using the -inf hack');
    nospike_term = -inf;
    %cifs blow up(usually because the k filters are too big)
end
log_prob = spike_term + nospike_term;


log_prob_full =  log_prob  + log(Trainpars.dt) * ( sum(Trainpars.logicalspike_microbin_Home) );
if (log_prob_full > 0)
    fprintf('Warning: encountered nonnegative log-likelihood value %f\n',log_prob_full);
end

    
end




