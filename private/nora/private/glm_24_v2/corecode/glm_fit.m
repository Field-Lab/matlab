function [p_opt,basepars,f_opt] = glm_fit(cid,basepars,spikes_lng,tstim,fn_save,gopts,X)
% This is a rip-off version of chaitu's mutiple_neuron.m.  Here the code is
% simplified to fit parameters of a *single* neuron (i.e., no coupling term
% between neurons).
% edoi@salk.edu, 2012-01-16,2012-04-16.

filt_mode = basepars.filtermode;
stimpars.dt = tstim; % this is the time (in sec) each stimulus frame is presented for.
basepars.stim_n = basepars.n; % no difference between stim_n and n.   x,y pts in stumulus space about the ROI

% Number of timebins used to store postspike filter
basepars.Mhist = floor((basepars.fac*basepars.Mk)*basepars.ps_fratio); % on fine timescale ~ 6/factor m(s)^-1
basepars.Nneurons = 1; % population size in GLM
trainpars.baseneuron_idx = 1; % index of the neuron(s) being trained simultaneously
trainpars.dt = stimpars.dt / basepars.fac; % time resolution of spikes

% linear filter basis -- not used.
basepars.ktime_basis = []; basepars.kspace_basis = [];

% ps basis parameters     AH  literally just parameters or basis
% construction
alpha_ps = basepars.alpha_ps; % starting point for PS
beta_ps  = (basepars.Mhist-1)*trainpars.dt; % ending point for ps filters

% History filters (postspike and coupling)
bstretch = basepars.bstretch;
if (basepars.nofilters_postspike > 0)
   basepars.postspike_basis = create_histbasis(alpha_ps,beta_ps,bstretch,...
      basepars.nofilters_postspike,basepars.Mhist,trainpars.dt,'L2',basepars.ps_spacing);
else
   basepars.postspike_basis = [];
end
% Coupling filters
basepars.coupling_basis = [];

% Obtain the spike times for each cell that is being used in training and
% store in a cell array as well as in a binary matrix
spikeT = basepars.fac*basepars.maxt;  % this is in (1/50)th units per frame

if (~isfield(basepars,'frame_offset'))
   stim_offset = 0;
else
   stim_offset = stimpars.dt*basepars.frame_offset;
end

trainpars.D = sparse(logical(false(spikeT,1)));  %in some sense initialize b4 puling spikes
duration = stimpars.dt * basepars.maxt; % amount of time in seconds of stimulus segment for training

%%% sp_times  is not really used,  but does give time in seconds for the
%%% conatenated blocks of interest
%%% trainpars.D  has the action, in the fine time scale logical 
[trainpars.sp_times,trainpars.D,trainpars.negSpikes] = pull_spikes(spikes_lng,cid,cid,spikeT,duration,trainpars.dt,stim_offset);

% Compute the postspike basis gradients (convolve each spike train with the basis fns)
fprintf('Computing the spike train convolutions...');
spikeT_offset = basepars.frame_offset*basepars.fac;
trainpars.psbasisGrad = grad_basis([trainpars.negSpikes; trainpars.D],basepars.postspike_basis,spikeT_offset);
trainpars.cpbasisGrad = []; % note: coupling is not examined here.
fprintf('done.\n');

%-- set the stimulus
frame_idx = basepars.frame_offset+1:basepars.frame_offset+basepars.maxt;    
stimpars.x = X(:,frame_idx); clear X

%%
p0 = basepars.p0;
basepars.crop_idx = (1:basepars.n)';
timeoffset = basepars.frame_offset*stimpars.dt; % offset time in seconds
fprintf('Using frames %d-%d spanning the time range %f-%f sec.\n',...
   basepars.frame_offset+1,basepars.frame_offset+basepars.maxt,timeoffset,timeoffset+duration);

%-- Fitting
switch filt_mode
   case 'rk2'
      [p_opt f_opt g_opt,H_opt,exitflag,output] = train_glm2(p0,basepars,stimpars,trainpars,gopts,1);
   case 'raw'  % raw = nonsep; pixel- and frame-based filter, i.e., no assumption on K.
      trainpars.lgrad = nonsep_lgrad(basepars,stimpars,trainpars);
      [p_opt,f_opt,g_opt,H_opt,exitflag,output] = train_glm_nonsep(p0,basepars,stimpars,trainpars,gopts);
end

%% save & print
opt_param.p  = p_opt;
opt_param.p0 = p0;
opt_param.f  = f_opt;
opt_param.g  = g_opt;
opt_param.H  = H_opt;
opt_param.exitflag = exitflag;
opt_param.output   = output;

stimpars = rmfield(stimpars,'x');
trainpars = rmfield(trainpars,'psbasisGrad');
if(isfield(trainpars,'lgrad'))
   trainpars = rmfield(trainpars,'lgrad');
end


%final outputs into the folder fit
save(fn_save,'basepars','opt_param','stimpars','trainpars','tstim','cid');
print_glm_fit(cid,opt_param,basepars,trainpars,fn_save);

