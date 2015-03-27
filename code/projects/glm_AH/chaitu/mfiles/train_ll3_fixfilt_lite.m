% Spike-driven version

% Function that computes the log likelihood of the data D for the parameters p.
% Returns the log probability of the data, the CIFS, and the filtered stimulus (under current parameter settings).
% The stimulus and the spike train can have different resolutions the dt should be specified correctly in pars and stimpars, resp.

function [lg_p input_term cifs kx lg_p_breakdown] = train_ll3_fixfilt_lite(p,basepars,stimpars,trainpars)

N = basepars.Nneurons;
Neff = size(trainpars.D,2);

T = size(stimpars.x,2);
numpars_perneuron = get_npars(basepars,Neff);

% Recover the base firing rates
boffset = get_pars_idx(basepars,1,Neff,'b');
b = p(boffset(1):numpars_perneuron:end)'; % 1  x N

% Recover stimulus filter norms
knormoffset = get_pars_idx(basepars,1,Neff,'k');
knorms = p(knormoffset(1):numpars_perneuron:end);
% edoi, editing, 2011-01-04.

kx = trainpars.kx * diag(knorms);
offset = get_pars_idx(basepars,1,Neff,'ps');
offset = offset(1)-1;

input_term = kx + repmat(b,T,1); % T x N

% Postspike filters corresponding to these parameters (dim: basepars.Mhist
% x neurons)

% Recover indices for the postpsike filters
add_term = repmat(numpars_perneuron,basepars.nofilters_postspike,N)*diag(0:(N-1)); % basepars.nofilters_postspike x N
idx = repmat((offset+1:offset+basepars.nofilters_postspike)',1,N) + add_term; % basepars.nofilters_postspike x N
ps_weights = reshape(p(idx(:)),basepars.nofilters_postspike,N);
PS = basepars.postspike_basis*ps_weights; % basepars.Mhist x N
offset = offset + basepars.nofilters_postspike;

% We need to "upsample" the input_term inorder to match the temporal
% resolution specified by the spike data (trainpars.dt)

factor = basepars.fac; %stimpars.dt / trainpars.dt;

% Just repeat every row of input term factor times
input_term = reprows(input_term,factor);

% Check that the lengths are matching up
if (factor < 1 || size(trainpars.D,1) ~= size(input_term,1))
    fprintf('ERROR: timebin sizes are still unmatched!!!!');
    return;
end

% Iterate through spikes instead of time

% Determine stimulus offset
if(~isfield(basepars,'frame_offset'))
    basepars.frame_offset = 0;
end

input_term = [zeros(basepars.frame_offset*basepars.fac,size(input_term,2)); input_term];
T2 = size(input_term,1);

if (~isfield(trainpars,'negSpikes'))
    trainpars.negSpikes = [];
end

spike_idx = find([trainpars.negSpikes; trainpars.D]);
D_idx = [(mod((spike_idx(:)-1),T2)+1) (floor((spike_idx(:)-1)/T2)+1)]; % spike list in special format - first column is timebin index, second col is neuron number

1;

% MEX fn to compute the log cifs
if (Neff > 1)
    input_term = ll_eval_mex(D_idx,input_term,PS,CP,trainpars.baseneuron_idx(:)); % time x N
else
    input_term = ll_eval_mex(D_idx,input_term,PS); % time x N
end

input_term = input_term(basepars.fac*basepars.frame_offset+1:end,:);


cifs = full(basepars.Nstep(input_term));
if (sum(cifs(:)<eps)>0)
    cifs = max(eps,cifs);
end

if (sum(cifs(:)>(1/trainpars.dt))>0)
    cifs = min(1/trainpars.dt,cifs);
end

% Determine relevant indices for training
rel_idx = get_train_idx(basepars,'spike');

spike_term =  sum(sum(log(cifs(rel_idx)).*trainpars.D(rel_idx,trainpars.baseneuron_idx)));

% Works only for exponential nonlinearity!
%spike_term = sum(sum(input_term.*trainpars.D(:,trainpars.baseneuron_idx)));

nospike_term = -trainpars.dt*sum(sum(cifs(rel_idx))); % approximation for trainpars.dt small


if (nargout > 3)
    spike_term_separate = sum(log(cifs(rel_idx)).*trainpars.D(rel_idx,trainpars.baseneuron_idx))'; % N x 1
    %spike_term_separate = sum(input_term.*trainpars.D(:,trainpars.baseneuron_idx))'; % N x 1
    nospike_term_separate = -trainpars.dt.*sum(cifs(rel_idx))'; % N x 1;
end

%nospike_term  = full(sum(sum(~trainpars.D .* log(1- trainpars.dt.*cifs))));

if (abs(max(cifs(:))) > 10^6)
   fprintf('Warning: large max cif value: logC = %f. spike_term=%f nospike_term=%f\n',full(log(max(cifs(:)))),full(spike_term),full(nospike_term));
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
lg_p = spike_term + nospike_term;
% We can add the constant term of the likelihood if necessary, just to be consistent - this
% is not needed for optimzation.
lg_p = lg_p + sum(sum(trainpars.D(rel_idx,trainpars.baseneuron_idx))).*log(trainpars.dt);
if (lg_p > 0)
    fprintf('Warning: encountered nonnegative log-likelihood value %f\n',lg_p);
    1;
end

if (nargout > 3)
   lg_p_breakdown  = spike_term_separate + nospike_term_separate;
   lg_p_breakdown = lg_p_breakdown + sum(trainpars.D(rel_idx,trainpars.baseneuron_idx))'.*log(trainpars.dt); % N x 1;
end
    
end