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
function [log_prob input_term cifs kx log_prob_breakdown] = train_ll3AH(p,Basepars,Stimpars,Trainpars,kxopt)

if isfield(Basepars,'Nneurons')
    N = Basepars.Nneurons;
else
    N = 1;
end


Neff = size(Trainpars.logicalspike_microbin_Home,2); 
Neighbors = length(Basepars.cp_Neighbors);


[nstim T] = size(Stimpars.movie_ROI);


numpars_perneuron = get_nparsAH(Basepars,Neighbors);
1;
%%% KX + MU    TERM OF THE GLM
% Recover the base firing rates
boffset = get_pars_idxAH(Basepars,1,Neff,'b');
b = p(boffset(1):numpars_perneuron:end)'; % 1  x N   %% looks like the mu term of the filter
%%% Run through stimulus filter portion  (no post spike or coupling yet)
if (~exist('kxopt','var'))
    if (isfield(Basepars,'XsqK') && Basepars.XsqK)
        %fprintf('Using square filter\n');
        1;
        [kx kxsq] = filterstimulus_trainAH(p,Basepars,Stimpars,Trainpars); % T x N
        kx = kx + kxsq;
    else
        kx = filterstimulus_trainAH(p,Basepars,Stimpars,Trainpars);
    end
else
    fprintf('train_ll3: Using the supplied kx of norm %f!\n',norm(kxopt));
    kx = kxopt;
end
1;
input_term = kx + repmat(b,T,1); % T x N   %%% AH i'm pretty sure this is analogous to the bas firing addition


% Check for extrinsic signal specification and add to input if necessary
if(isfield(Basepars,'extrinsic') && ~isempty(Basepars.extrinsic))
    input_term = input_term + Basepars.extrinsic(Basepars.frame_offset+1:Basepars.frame_offset+T,:);
end
if (isfield(Basepars,'ext_timepts') && ~isempty(Basepars.ext_timepts) > 0)
    
    1;
    % Assume that it's a band-limited function with
    % length(Basepars.ext_timepts) being the critical sampling - do
    % bandlimited-interpolation on these samples.
    input_term = input_term + geninterp(Basepars.interpfilt,Basepars.maxt,Basepars.ext_timepts,p(get_pars_idx(Basepars,1,Neff,'ext')));
    %sincinterp(makeaxis(Stimpars.dt,Basepars.maxt),Basepars.ext_timepts,p(get_pars_idx(Basepars,1,Neff,'ext')));
end


%%% POST SPIKE FILTERING %%%%

% Recover indices for the postpsike filters     %%% AH prolly makes more
% sense when we are using multiple neurons
if (Basepars.ps_filternumber > 0)
    psidx = get_pars_idxAH(Basepars,1,Neff,'ps');  %indexes the ps indices of p
    offset = psidx(1)-1;
    add_term = repmat(numpars_perneuron,Basepars.ps_filternumber,N)*diag(0:(N-1)); % Basepars.ps_filternumber x N
    idx = repmat((offset+1:offset+Basepars.ps_filternumber)',1,N) + add_term; % Basepars.ps_filternumber x N
    ps_weights = reshape(p(idx(:)),Basepars.ps_filternumber,N);
    PS = Basepars.ps_basis*ps_weights; % Basepars.Mhist x N    
else
    PS = zeros(Basepars.ps_timebins,1);
end

%%% COUPLING HAPPENS HERE !!!  

if (Neighbors >= 1 && Basepars.cp_filternumber > 0)
    cpidx = get_pars_idxAH(Basepars,1,Neff,'cp');
    offset = cpidx(1)-1;

    % Coupling filters corresponding to these parameters (dim: Basepars.Mhist x neuron1*neuron2)
    add_term = repmat(numpars_perneuron,(Neff-1)*Basepars.cp_filternumber,N)*diag(0:N-1); % (N-1)*nofitlers_coupling x N
    idx = repmat((offset+1:offset+(Neff-1)*Basepars.cp_filternumber)',1,N) + add_term; % (N-1)*Basepars.cp_filternumber x N
    c_weights = reshape(p(idx(:)),Basepars.cp_filternumber,N*(Neff-1)); % Basepars.cp_filternumber x N*(Neff-1)
    1;
    CP = Basepars.coupling_basis*c_weights;
    ncoupling = size(c_weights,1); 
end


 
%%%%%   UPSAMPLE THE INPUT_TERM TO MATCH TEMPORAL RES OF SPIKE DATA  %%%%
factor = Basepars.spikebins_perstimframe; %Stimpars.dt / Trainpars.dt;
% Just repeat every row of input term factor times
input_term = reprows(input_term,factor);  %% NOW ON MICROBIN TIME SCALE
%    DETERMINE IF THERE SHOULD BE ANY OFFSET
if(~isfield(Basepars,'frame_offset'))
    Basepars.frame_offset = 0;
end
if (isfield(Basepars,'analog_frame_offset') && ~isempty(Basepars.analog_frame_offset))
    numbinsadded = double(int32(Basepars.analog_frame_offset/Trainpars.dt));
else
    numbinsadded = Basepars.frame_offset*Basepars.spikebins_perstimframe;    
end
input_term = [zeros(numbinsadded,size(input_term,2)); input_term];
T2 = size(input_term,1);   %% NUMBER OF MICROBINS
% T EQUALS THE NUMBER OF STIM FRAMES
if (~isfield(Trainpars,'negSpikes'))
    Trainpars.negSpikes = [];
end


%%%  THIS IS THE INDEX TIME OF OUR SPIKES %%%
%%%  IMPORTANT INPUT INTO  LL_EVAL_MEX    %%%
spike_idx = find([Trainpars.negSpikes_Home; Trainpars.logicalspike_microbin_Home]);
D_idx = [(mod((spike_idx(:)-1),T2)+1) (floor((spike_idx(:)-1)/T2)+1)]; % spike list in special format - first column is timebin index, second col is neuron number


%% NEED TO THINK ABOUT STRUCTURE WHEN WE ADD COUPLING 
if (size(input_term,1) ~= (size(Trainpars.negSpikes,1)+size(Trainpars.logicalspike_microbin_Home,1)))
    error('input term has length %d and does not match concatenated spike train length %d',size(input_term,1),(size(Trainpars.negSpikes,1)+size(Trainpars.logicalspike_microbin_Home,1)));
end



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MEX fn to compute the log cifs %%%%%%%%%%%%%%%%
%%%%% Matlab code below is almost as fast %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% take less than half a second

if (Neighbors > 1)
    input_term = ll_eval_mex(D_idx,input_term,PS,CP,Trainpars.baseneuron_idx(:)); % time x N
    
else
    input_term = ll_eval_mex(D_idx,input_term,PS); % time x N
end
input_term = input_term(numbinsadded+1:end,:);

% -------------------------------------------------
if 0 % this is the MATLAB version - much slower!   %% NOT MUCH SLOWER... MIGHT AS WELL USE THIS IF NOT COUPLED
    tic
    for i=1:length(spike_idx)
        1;
        % Recover the time and the neuron corresponding to this spike
        n2 = floor((spike_idx(i)-1)/T2)+1;
        t = mod((spike_idx(i)-1),T2)+1;
        %[t n2] = ind2sub(size(Trainpars.logicalspike_microbin_Home),spike_idx(i)); %remember, n2 is the neuron affecting another neuron n1
        
        timerange = (t+1):min(t+Basepars.ps_timebins,T2); % Explicit time range that is affected
        numtoupdate = length(timerange); % Number of future timebins this spike will affect
        
        % Compute the recurrent (postspike effects) contribution of this spike
        input_term(timerange,n2) = input_term(timerange,n2) + PS(1:numtoupdate,n2);
        
        % Network contributions (coupling effects) of this spike
        if (N > 1)
            for n1=[1:n2-1 n2+1:N] % n1 is the neuron BEING AFFECTED - we should update its input term
                n2rev = n2;
                if (n2 > n1)
                    n2rev = n2rev - 1;
                end
                colnum = (N-1)*(n1-1)+(n2rev-1)+1;
                %fprintf('trainll2: Checking the effect of neuron %d on neuron %d by looking at column %d\n',n2,n1,colnum);
                update_term = C(1:numtoupdate,colnum);
                input_term(timerange,n1) = input_term(timerange,n1) + update_term;
                
            end
        end
    end
    toc
end
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