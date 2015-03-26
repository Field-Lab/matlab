%%% THis one is just designed to be flexible to include coupling as well
%% Default is matlab based code, not  C based code 
%%% WORKS AS OF 10/10/2012 

% Spike-driven version
%%%%%%%%% 
%%%%%%%%% CALLS ll_eval_mex %%%%%%%%%%%%
%% BUT DOENS'T HAVE TO.. MATLABECODE IS JUST TWICE SLOWER %%%

%INPUT_TERM MEANS LCIF FOR NOW

% FUNCTINO CALLS: FILTERSTIMULUSTRAIN, GET_NPARS, GET_PARS_IDX
% ALSO CALLS GENINTERP (BUT NOT IN MY PARAM RANGE)


% Returns the log probability of the data, the CIFS, and the filtered
% stimulus (under current parameter settings).

% The stimulus and the spike train can have different resolutions
% the dt should be specified correctly in pars and Stimpars, resp.

%% GENERATE  A KX A PS AND A CP  AND COMPUTE FROM THERE ?? 


%% NEED TO KEEP TRACK OF INPUT_TERM .. BECOMES THE FINAL LCIFS
function [log_prob input_term cifs kx log_prob_breakdown] = train_ll3_cpAH(p,Basepars,Stimpars,Trainpars,kxopt)

if ~isfield(Basepars, 'LL_eval_C')
    Basepars.LL_eval_C = false 
end

nNeighbors = length(Basepars.cp_Neighbors);
[spacePixels stimFrames] = size(Stimpars.movie_ROI);
numpars = length(Basepars.p0);
microbins = stimFrames * Basepars.spikebins_perstimframe;


%%% KX + MU    TERM OF THE GLM
mu_index = get_pars_idx2AH(Basepars,'mu');
mu = p(mu_index);% 1  x N   %% looks like the mu term of the filter
%%% Run through stimulus filter portion  (no post spike or coupling yet)
if (~exist('kxopt','var'))
    if (isfield(Basepars,'XsqK') && Basepars.XsqK)
        %fprintf('Using square filter\n');
        1;
        [kx kxsq] = filterstimulus_trainAH(p,Basepars,Stimpars); % T x N
        kx = kx + kxsq;
    else
        kx = filterstimulus_trainAH(p,Basepars,Stimpars);
    end
else
    fprintf('train_ll3: Using the supplied kx of norm %f!\n',norm(kxopt));
    kx = kxopt;
end
input_term = kx + repmat(mu,stimFrames,1); % T x N   %%% AH i'm pretty sure this is analogous to the bas firing addition


%%% POST SPIKE FILTERING %%%%

if (Basepars.ps_filternumber > 0)
    psidx = get_pars_idx2AH(Basepars,'ps');  %indexes the ps indices of p
    ps_weights = p(psidx(:));
    PS = Basepars.ps_basis*ps_weights; % Basepars.Mhist x N    
else
    PS = zeros(Basepars.ps_timebins,1);
end
%PSupdates = (ps_weights') *Trainpars.psbasisGrad; 

%%% COUPLING HAPPENS HERE !!!  %%%%%

if (Basepars.Coupling && nNeighbors> 0)
    cpidx = get_pars_idx2AH(Basepars,'cp');
    cp_weights = reshape(p(cpidx(:)),Basepars.cp_filternumber,nNeighbors); % THIS NEEDS TO BE CONSISTENT IN WHAT P MAPS TO!! %%%
    CP = Basepars.cp_basis*cp_weights;  %% CP is a matrix !!
end

 
%%%%% PUT INPUT_TERM (KX+MU) ONTO MICROBIN TIMESCALE  %%%%
input_term = reprows(input_term,Basepars.spikebins_perstimframe); 

%%%%    DETERMINE IF THERE SHOULD BE ANY TIMING OFFSET  %%%%
if(~isfield(Basepars,'frame_offset'))
    Basepars.frame_offset = 0;
end
if (isfield(Basepars,'analog_frame_offset') && ~isempty(Basepars.analog_frame_offset))
    numbinsadded = double(int32(Basepars.analog_frame_offset/Trainpars.dt));
else
    numbinsadded = Basepars.frame_offset*Basepars.spikebins_perstimframe;    
end
input_term = [zeros(numbinsadded,size(input_term,2)); input_term];


if (~isfield(Trainpars,'negSpikes_Home'))
    Trainpars.negSpikes_Home = [];
end
if (~isfield(Trainpars,'negSpikes_Neighbors'))
    Trainpars.negSpikes_Neighbors = [];
end

%%%  THIS IS THE INDEX TIME OF OUR SPIKES %%%
%%%  IMPORTANT INPUT INTO  LL_EVAL_MEX    %%%
spike_idx_Home = find([Trainpars.negSpikes_Home; Trainpars.logicalspike_microbin_Home]);
D_idx_Home = [(mod((spike_idx_Home(:)-1),microbins)+1) (floor((spike_idx_Home(:)-1)/microbins)+1)]; 
% spike list in special format - first column is timebin index, second col
% is neuron number
%%%  IMPORTANT INPUT INTO  LL_EVAL_MEX    %%%
spike_idx_Neighbors = find([Trainpars.negSpikes_Neighbors; Trainpars.logicalspike_microbin_Neighbors]);
D_idx_Neighbors = [(mod((spike_idx_Neighbors(:)-1),microbins)+1) (floor((spike_idx_Neighbors(:)-1)/microbins)+1)]; 

if (size(input_term,1) ~= (size(Trainpars.negSpikes_Home,1)+size(Trainpars.logicalspike_microbin_Home,1)))
    error('input term has length %d and does not match concatenated spike train length %d',size(input_term,1),(size(Trainpars.negSpikes,1)+size(Trainpars.logicalspike_microbin_Home,1)));
end
%% NEED TO THINK ABOUT STRUCTURE WHEN WE ADD COUPLING 


%%%% Modify mu + kx  by PS and CP%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MEX fn to compute the log cifs %%%%%%%%%%%%%%%%
%%%%% Matlab code below is almost as fast %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% take less than half a second
%{
if Basepars.LL_eval_C
    if (nNeighbors > 1)
        input_term = ll_eval_mex(D_idx,input_term,PS,CP,Trainpars.baseneuron_idx(:)); % time x N

    else
        input_term = ll_eval_mex(D_idx_Home,input_term,PS); % time x N
    end
    input_term = input_term(numbinsadded+1:end,:);

% -------------------------------------------------
else  % this is the MATLAB version - much slower!   %% NOT MUCH SLOWER... MIGHT AS WELL USE THIS 
        
    %%%% PS Component %%%%%%
    PSupdate = zeros(microbins,1);
    for i=1:length(spike_idx_Home)
        t = spike_idx_Home(i);
        timerange = (t+1):min(t + Basepars.ps_timebins,microbins); % Explicit time range that is affected
        numtoupdate = length(timerange); % Number of future timebins this spike will affect      
        % Compute the recurrent (postspike effects) contribution of this spike
        PSupdate(timerange) = PSupdate(timerange) + PS(1:numtoupdate);
    end 
    input_term = input_term + PSupdate;
    
    %%%% CP Component %%%%%
	if (Basepars.Coupling && nNeighbors> 0)
        CPupdates = zeros(microbins, nNeighbors);
        for i = 1:size(D_idx_Neighbors,1)
            t    = D_idx_Neighbors(i,1);
            cell = D_idx_Neighbors(i,2);
            timerange = (t+1):min(t + Basepars.cp_timebins,microbins); 
            numtoupdate = length(timerange);
            
            CPupdates(timerange,cell) = CPupdates(timerange,cell) + CP(1:numtoupdate,cell);
        end
        CPupdate = sum(CPupdates,2);
        input_term = input_term + CPupdate;
    end
end
%}
% --------------------------------------------

%%% AT THIS POINT THE INPUT TERM IS NOW EQUAL TO THE LCIFS
cifs = full(Basepars.Nstep(input_term));
if (sum(cifs(:)<eps)>0)
    cifs = max(eps,cifs);
end

if (sum(cifs(:)>(1/Trainpars.dt))>0)
    cifs = min(1/Trainpars.dt,cifs);
end

% Determine relevant indices for training    get_train_idx   is short and
% will just implement it here
%%rel_idx = get_train_idx(Basepars,'spike');
res = 'spike';
    if (~isfield(Basepars,'leaveout_idx'))
        idx = 1:Basepars.maxt;
    else
        idx = setdiff(1:Basepars.maxt,Basepars.leaveout_idx);
    end
    idx = idx(:);

    if  strcmp(res,'spike')
        idx = reshape(repmat(Basepars.spikebins_perstimframe*(idx'-1),Basepars.spikebins_perstimframe,1)...
            + repmat((1:Basepars.spikebins_perstimframe)',1,length(idx)),length(idx)*Basepars.spikebins_perstimframe,1);
    end
rel_idx  = idx;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   COMPUTING THE LOG PROB   %%%%%%%%%%%%%%
spike_term =  sum(sum(log(cifs(rel_idx)).*Trainpars.logicalspike_microbin_Home(rel_idx,Trainpars.baseneuron_idx)));
% Works only for exponential nonlinearity!
%spike_term = sum(sum(input_term.*Trainpars.logicalspike_microbin_Home(:,Trainpars.baseneuron_idx)));
nospike_term = -Trainpars.dt*sum(sum(cifs(rel_idx))); % approximation for Trainpars.dt small

%%% STUFF THAT DOESN'T SEEM RELEVANT YET   %%%%%%%%%
if (nargout > 3)
    spike_term_separate = sum(log(cifs(rel_idx)).*Trainpars.logicalspike_microbin_Home(rel_idx,Trainpars.baseneuron_idx))'; % N x 1
    %spike_term_separate = sum(input_term.*Trainpars.logicalspike_microbin_Home(:,Trainpars.baseneuron_idx))'; % N x 1
    nospike_term_separate = -Trainpars.dt.*sum(cifs(rel_idx))'; % N x 1;
end
%nospike_term  = full(sum(sum(~Trainpars.logicalspike_microbin_Home .* log(1- Trainpars.dt.*cifs))));
if (abs(max(cifs(:))) > 10^6)
    fprintf('Warning: large max cif value: logC = %f. spike_term=%f nospike_term=%f\n',full(log(max(cifs(:)))),full(spike_term),full(nospike_term));
    %fprintf('b=%f,ks1=%f,ks2=%f,kt1=%f,kt2=%f,ps=%f\n',p(1),norm(p(1+1:1+Basepars.n)),norm(p(1+Basepars.n+Basepars.nofilters_k+1:1+Basepars.n+Basepars.nofilters_k+Basepars.n)),...
    %                                                    norm(p(1+Basepars.n+1:1+Basepars.n+Basepars.nofilters_k)),norm(p(1+2*Basepars.n+Basepars.nofilters_k+1:1+2*(Basepars.n+Basepars.nofilters_k))),...
    %                                                    norm(p(end-Basepars.ps_filternumber+1:end)));
end
% This is a hack - if the cifs are blowing up during optimization (usually because the k
% filters are too big), we penalize by -inf
if (~isreal(nospike_term))
    fprintf('Warning: No spike is imaginary - using the -inf hack');
    nospike_term = -inf;
end
if (nargout > 3)
    imag_idx = find(abs(imag(nospike_term_separate)) > eps);
    nospike_term_separate(imag_idx) = -inf;
end
% Determine the idx to be used for training (there may be a leaveout_idx specified)
log_prob = spike_term + nospike_term;
% We can add the constant term of the likelihood if necessary, just to be consistent - this
% is not needed for optimzation.
log_prob = log_prob + sum(sum(Trainpars.logicalspike_microbin_Home(rel_idx,Trainpars.baseneuron_idx))).*log(Trainpars.dt);
if (log_prob > 0)
    fprintf('Warning: encountered nonnegative log-likelihood value %f\n',log_prob);
    1;
end

if (nargout > 3)
   log_prob_breakdown  = spike_term_separate + nospike_term_separate;
   log_prob_breakdown = log_prob_breakdown + sum(Trainpars.logicalspike_microbin_Home(rel_idx,Trainpars.baseneuron_idx))'.*log(Trainpars.dt); % N x 1;
end
    
end






%{
% Check for extrinsic signal specification and add to input if necessary
if(isfield(Basepars,'extrinsic') && ~isempty(Basepars.extrinsic))
    input_term = input_term + Basepars.extrinsic(Basepars.frame_offset+1:Basepars.frame_offset+T,:);
end
if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts) > 0)
    % Assume that it's a band-limited function with
    % length(Basepars.ext_timepts) being the critical sampling - do
    % bandlimited-interpolation on these samples.
    input_term = input_term + geninterp(Basepars.interpfilt,Basepars.maxt,Basepars.ext_timepts,p(get_pars_idx(Basepars,1,Neff,'ext')));
    %sincinterp(makeaxis(Stimpars.dt,Basepars.maxt),Basepars.ext_timepts,p(get_pars_idx(Basepars,1,Neff,'ext')));
end
}
