function [spikes cif hist_term xout kxout] = simulate_trains2(pars,opts,stim,Ntrials,chunk_size)

if ~exist('chunk_size','var')
    chunk_size = -1;
end

% For loop instead of 3D matrix

spikes = ones(pars.maxt,Ntrials,pars.Nneurons) > 2;
cif = zeros(pars.maxt,pars.Nneurons);
hist_term = zeros(pars.maxt,Ntrials,pars.Nneurons);

% Filter the stimulus - kx will have dims (time x neurons)
if (isnumeric(stim))
    if (opts.filterstimulus)
        kx = filterstimulus(pars,opts,stim);
    else % Turn filtering off - just take mean in spatial dimension and replicate for each neuron.
        kx = repmat(mean(stim,1)',1,pars.Nneurons);
    end
end

if 0%opts.plotfilteredstimuli
    figure, plot(pars.t,kx);
    xlabel('Time');
    ylabel('Filtered stimulus');
end

if (nargout > 3)
    xout = zeros(pars.n,pars.maxt,Ntrials);
end

if (nargout > 4)
    kxout = zeros(pars.maxt,pars.Nneurons,Ntrials);
end

for j=1:Ntrials
    if (mod(j,100) == 0)
        fprintf('Completed %d trials.\n',j);
    end
   
    
    if(isstruct(stim))
        % Generate a GWN stimulus
        xj = stim.sig.*randn(pars.n,pars.maxt)+stim.mu;
        if (opts.filterstimulus)
            kx = filterstimulus(pars,opts,xj);
        else
            kx = repmat(mean(xj,1)',1,pars.Nneurons);
        end
        
        if (nargout > 3) %store the stimulus as well
            xout(:,:,j) = xj;
        end
        
        if (nargout > 4) % store the filtered stimulus as well
            kxout(:,:,j) = kx;
        end
        
    end
    [spikej inputj hist_termj] = simulate_trains_fast(pars,opts,kx,chunk_size);
    spikes(:,j,:) = spikej;
    cif = cif + (1/Ntrials).*pars.Nstep(inputj);
    hist_term(:,j,:) = hist_termj;
end
