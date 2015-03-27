function [spikes lcifs hist_term] = simulate_reptrials(basepars,X,ntrials,b,K,PS,dt)

if (~isnumeric(X))
    if (strcmp(X,'bwn'))
        opt_pars.mode = 'switch';
        opt_pars.mu1 = 0.5;
        opt_pars.mu2 = 0.5;
        opt_pars.sig1 = 0.0;
        opt_pars.sig2 = 0.48;
        opt_pars.noisetype = 'bwn';
        opt_pars.period = 5;
        basepars.maxt = 2000;%ceil(opt_pars.period*10/(dt*basepars.fac));
        basepars.crop_idx{1} = (1:basepars.n)';
    end
else
    opt_pars.X = X;
    basepars.maxt = size(X,2);
end
for k=1:ntrials
    [spikes(:,k) lcifs(:,k) hist_term(:,k)] = simulate_model_nopars(basepars,b,K,PS,[],[],[],dt,1,opt_pars);
end
    