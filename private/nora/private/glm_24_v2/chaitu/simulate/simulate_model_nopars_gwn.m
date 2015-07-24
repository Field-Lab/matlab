function [spikes lcifs hist_term kx X] = simulate_model_nopars_gwn(base_params,b,K,PS,spike_dt,mu,p,ntrials,X)
    spikes = [];
    retries = 1;
    while (isempty(spikes))
        if (retries > 1)
            1;
            fprintf('Failed to simulate for this GWN stimulus. Resampling stimulus attempt number %d\n',retries);
            RandStream.setDefaultStream(RandStream.create('mt19937ar','Seed',cputime));
        end
        if (nargin < 9)
            X = create_gwn(mu,p,base_params.stim_dt,size(K,2),base_params.stimT,[0 inf]);
        else
            fprintf('simulate_model_nopars_gwn: Using supplied stimulus X...\n');
        end
        base_params.maxt = size(X,2);    
        1;
        tic; [spikes lcifs hist_term kx] = simulate_model_nopars(base_params,b,K',PS,[],X,[],spike_dt,ntrials); toc;
        retries = retries+1;
    end

    
    1;