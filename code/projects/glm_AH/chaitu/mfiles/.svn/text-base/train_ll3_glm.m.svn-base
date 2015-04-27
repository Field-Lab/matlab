% Returns the log probability of the data, the CIFS, and the filtered stimulus
% (under current parameter settings).

% The stimulus and the spike train can have different resolutions.
% The "dt" should be specified correctly in pars and stimpars, resp.

% this code is based on train_ll3_lite.m that is specialized for vGLM.

function [lg_p,input_term,cifs,az,lg_p_breakdown] = train_ll3_glm(p,basepars,stimpars,trainpars)
% Function that computes the log likelihood of the data D for the parameters p

% "input_term" aka "lcifs" in the other functions.

N = 1; % = 1 in lite
Neff = 1; % = 1 in lite
basepars.Nneurons = 1;

[~, T] = size(stimpars.z);
numpars_perneuron = get_npars(basepars,Neff);

%-- baseline
% boffset = get_pars_idx(basepars,1,Neff,'b');
% b = p(boffset(1):numpars_perneuron:end)'; % 1  x N
b = p(1);

% Recover stimulus filter norms
%knormoffset = get_pars_idx(basepars,1,Neff,'k');
%knorms = p(knormoffset(1):numpars_perneuron:end);

% edoi, editing, 2011-01-04.
% kx = trainpars.kx * diag(knorms);
% offset = get_pars_idx(basepars,1,Neff,'ps');
% offset = offset(1)-1;
%kidx = 2:(basepars.stim_n+1);

%keyboard % 2012-02-15

az = (stimpars.z)';

% %-- filtered stimulus
% if (~exist('kxopt','var'))
%    % Filter the stimulus - use the HST constant!
%    kx = filterstimulus_train(p,basepars,stimpars,trainpars);
%    
% else
%    fprintf('train_ll3: Using the supplied kx of norm %f!\n',norm(kxopt));
%    kx = kxopt;
%    
% end

%-- sum of the above two
input_term = az + repmat(b,T,1); % T x N


%-- Recover indices for the postpsike filters
if (basepars.nofilters_postspike > 0)
   %psidx = get_pars_idx(basepars,1,Neff,'ps');
   psidx = basepars.stim_n+2:length(p); % change later
   offset = psidx(1)-1;
   add_term = repmat(numpars_perneuron,basepars.nofilters_postspike,N)*diag(0:(N-1)); % basepars.nofilters_postspike x N
   idx = repmat((offset+1:offset+basepars.nofilters_postspike)',1,N) + add_term; % basepars.nofilters_postspike x N
   ps_weights = reshape(p(idx(:)),basepars.nofilters_postspike,N);
   PS = basepars.postspike_basis*ps_weights; % basepars.Mhist x N
else
   PS = zeros(basepars.Mhist,1);
end

% note: assume no coupling in lite version (edoi)

% We need to "upsample" the input_term inorder to match the temporal resolution specified by the spike data (trainpars.dt)
factor = basepars.fac; 
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
if (isfield(basepars,'analog_frame_offset') && ~isempty(basepars.analog_frame_offset))
    numbinsadded = double(int32(basepars.analog_frame_offset/trainpars.dt));
else
    numbinsadded = basepars.frame_offset*basepars.fac;
end
input_term = [zeros(numbinsadded,size(input_term,2)); input_term];

T2 = size(input_term,1);

if (~isfield(trainpars,'negSpikes'))
    trainpars.negSpikes = [];
end

spike_idx = find([trainpars.negSpikes; trainpars.D]);
D_idx = [(mod((spike_idx(:)-1),T2)+1) (floor((spike_idx(:)-1)/T2)+1)];
% spike list in special format - first column is timebin index, second col is neuron number

if (size(input_term,1) ~= (size(trainpars.negSpikes,1)+size(trainpars.D,1)))
   error('input term has length %d and does not match concatenated spike train length %d',...
      size(input_term,1),(size(trainpars.negSpikes,1)+size(trainpars.D,1)));
end

%-- MEX fn to compute the (log) cifs % see train_ll3.m for matlab implementation
input_term = ll_eval_mex(D_idx,input_term,PS); % time x N
input_term = input_term(numbinsadded+1:end,:); % edoi: final input term
cifs = full(basepars.Nstep(input_term));
if (sum(cifs(:)<eps)>0)
   cifs = max(eps,cifs);
end

if (sum(cifs(:)>(1/trainpars.dt))>0)
   cifs = min(1/trainpars.dt,cifs);
end

%-- Determine relevant indices for training
rel_idx = get_train_idx(basepars,'spike');
%spike_term = sum(sum(log(cifs(rel_idx)).*trainpars.D(rel_idx,trainpars.baseneuron_idx)));
spike_term = sum(sum(log(cifs(rel_idx)).*trainpars.D(rel_idx)));

nospike_term = -trainpars.dt*sum(sum(cifs(rel_idx))); % approximation for trainpars.dt small

if (nargout > 3)
   %spike_term_separate = sum(log(cifs(rel_idx)).*trainpars.D(rel_idx,trainpars.baseneuron_idx))'; % N x 1
   spike_term_separate = sum(log(cifs(rel_idx)).*trainpars.D(rel_idx))'; % N x 1
   nospike_term_separate = -trainpars.dt.*sum(cifs(rel_idx))'; % N x 1;
end

if (abs(max(cifs(:))) > 10^6)
   fprintf('Warning: large max cif value: logC = %f. spike_term=%f nospike_term=%f\n',...
      full(log(max(cifs(:)))),full(spike_term),full(nospike_term));
end

% This is a hack - if the cifs are blowing up during optimization (usually because the k filters are too big),
% we penalize by -inf
if (~isreal(nospike_term))
   fprintf('Warning: No spike is imaginary - using the -inf hack');
   nospike_term = -inf;
end
if (nargout > 3)
   imag_idx = find(abs(imag(nospike_term_separate)) > eps);
   nospike_term_separate(imag_idx) = -inf;
end
        
%-- Determine the idx to be used for training (there may be a leaveout_idx specified)
lg_p = spike_term + nospike_term;
% We can add the constant term of the likelihood if necessary, just to be consistent - this
% is not needed for optimzation.
%lg_p = lg_p + sum(sum(trainpars.D(rel_idx,trainpars.baseneuron_idx))).*log(trainpars.dt);
lg_p = lg_p + sum(sum(trainpars.D(rel_idx))).*log(trainpars.dt);
if (lg_p > 0)
   fprintf('Warning: encountered nonnegative log-likelihood value %f\n',lg_p);
end

if (nargout > 3)
   lg_p_breakdown  = spike_term_separate + nospike_term_separate;
   %lg_p_breakdown = lg_p_breakdown + sum(trainpars.D(rel_idx,trainpars.baseneuron_idx))'.*log(trainpars.dt); % N x 1;
   lg_p_breakdown = lg_p_breakdown + sum(trainpars.D(rel_idx))'.*log(trainpars.dt); % N x 1;
end
    
end