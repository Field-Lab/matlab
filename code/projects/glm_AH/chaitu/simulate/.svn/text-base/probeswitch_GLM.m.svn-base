% Function simulates a luminance/contrast switching experiment on a GLM
% neuron

% Arguments:

% basepars - base parameters of GLM
% b,K,PS,CP - LIST model parameters (base rate, stim filt, postspike filt, coupling filt)
% dt,ntrials - timestep, # of trials
% opt_pars - stimulus switching information (pre and postswitch luminance and contrast)
% periods - different switching periods to try

% Returns:

% results - a matrix of structures indexed (stimulus pars) x (model pars)
%           containing the firing rates and spikes of the neurons



function results = probeswitch_GLM(basepars,b,K,PS,CP,dt,ntrials,opt_pars,periods)

if(~iscell(b))
    b = num2cell(b);
end

nsets = length(b);
nperiods = length(periods);
results = cell(nperiods,nsets);


for i=1:nperiods
    for j=1:nsets
        
        % Set up the number of frames to account for the period size
        if (isfield(opt_pars,'X'))
            nframes_perperiod = size(opt_pars.X,2);
        else
            nframes_perperiod = ceil(periods(i)/(dt*basepars.fac));
        end
        basepars.maxt = ceil(4.1*nframes_perperiod);
        basepars.n = size(K{j},1);
        basepars.crop_idx = (1:basepars.n)';
        opt_pars.period = nframes_perperiod;
        
        if (isempty(CP))
            CPj = [];
        else
            CPj = CP{j};
        end
        fprintf('Simulating (%d/%d,%d/%d) GLM with stimulus period %f and PS norm %f.\n',i,nperiods,j,nsets,periods(i),norm(PS{j}));
        1;
        [spikes lcifs] = simulate_model_nopars(basepars,b{j},K{j},PS{j},CPj,[],[],dt,ntrials,opt_pars);
        lcifs = squeeze(lcifs);
        
        % Take the relevant segment and average rates (disregard the 1st one)
        T = nframes_perperiod*basepars.fac;
        eff_end = T*floor(basepars.maxt/nframes_perperiod);
        idx = T+1:eff_end;
        %stim_idx = nframes_perperiod+1:eff_end/basepars.fac;
        ntrials_eff = (floor(basepars.maxt/nframes_perperiod)-1)*ntrials;

        1;
        
        results{i,j}.rates = exp(reshape(squeeze(lcifs(idx,:,:)),T,ntrials_eff));
        results{i,j}.spikes = reshape(squeeze(spikes(idx,:,:)),T,ntrials_eff);
        %results{i,j}.kx = reshape(squeeze(kx(stim_idx,:,:)),nframes_perperiod,ntrials_eff);
        %figure, plot(taxis,log(mean(exp(squeeze(lcifs)),2)));
        %hold on, plot(taxis,reprows(kx,basepars.fac),'r-')
    end
    
end

