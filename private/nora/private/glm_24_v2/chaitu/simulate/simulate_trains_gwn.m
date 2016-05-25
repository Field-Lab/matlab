function [spikes cif_stats hist_term] = simulate_trains_gwn(pars,opts,Ntrials,stimpars)

if ~exist('chunk_size','var')
    chunk_size = -1;
end

% For loop instead of 3D matrix

spikes = ones(pars.maxt,Ntrials,pars.Nneurons) > 2;
cif_stats.avg_cifs = zeros(pars.maxt,pars.Nneurons);
cif_stats.std_cifs = zeros(pars.maxt,pars.Nneurons);
hist_term = zeros(pars.maxt,Ntrials,pars.Nneurons);

for j=1:Ntrials
    if (mod(j,100) == 0)
        fprintf('Completed %d trials.\n',j);
    end

    % Create random stimulus from stimpars on each trial
    x = randn(pars.n,pars.maxt);
    x(:,1:stimpars.jump_idx) = x(:,1:stimpars.jump_idx).*stimpars.sig1 + stimpars.mu1;
    x(:,stimpars.jump_idx+1:end) = x(:,stimpars.jump_idx+1:end).*stimpars.sig2 + stimpars.mu2;
    % ---
    kx = filterstimulus(pars,opts,x);
    [spikej inputj hist_termj] = simulate_trains_fast(pars,opts,kx,chunk_size);
    spikes(:,j,:) = spikej;
    cif_stats.avg_cifs = cif_stats.avg_cifs + (1/Ntrials).*pars.Nstep(inputj);
    hist_term(:,j,:) = hist_termj;
end
