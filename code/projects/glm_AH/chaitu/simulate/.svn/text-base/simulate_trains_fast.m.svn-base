function [spikes input hist_term] = simulate_trains_fast(pars,opts,kx,chunk_size)

% Function for simulating spike trains given a filtered stimulus
% Args: parameters, options, stimulus x, number of trial Ntrials
% Returns:
% spikes - binary matrix (time x trial x neuron) of spikes 
% cif_stats - contains the average CIF and the standard deviation of the CIF for each time bin
% hist_term - h currents (time x trial x neuron) (for debugging mostly)



% Initialize the simulation variables
spikes = (zeros(pars.maxt,pars.Nneurons) > 1); % the spike trains (logical matrix)
input = repmat(pars.b',pars.maxt,1) + kx; % the input current (log of the CIF) pars.maxt x pars.Nneurons
hist_term = zeros(pars.maxt,pars.Nneurons); % the history terms
%DEBUGGING
%hist_test = zeros(pars.maxt,Ntrials,pars.Nneurons);




% Cap if necessary (only do this for the first time step - there is an
% error if we have to do this in the middle)
cap_idx = input(1,:) > log(1/pars.dt);
input(1,cap_idx) = log(1/pars.dt);
%if (sum(sum(sum(cap_idx))) > 0)
%    fprintf('Capped the initial CIF value for %d neurons to %d\n\n',sum(sum(sum(cap_idx))),1/pars.dt);
%end

% Time window to stimulate per iteration
DEFAULT_CHUNK_SIZE = floor(max(100,200/(exp(pars.b+mean(kx)))));

if (~exist('chunk_size','var') || chunk_size < 0)
    chunk_size = DEFAULT_CHUNK_SIZE;
end
%chunk_size
t = 0; % time index (not ACTUAL time!)

% Run simulation till end of time interval
while (t < pars.maxt)
    timerange = (t+1):min(t+chunk_size,pars.maxt);
    
    % Generate a random Exp(1) value
    isi = exprnd(1,1,pars.Nneurons);
    % Find the first index at which the integrated cifs exceed the isi    
    cif_chunk = cumsum(pars.dt .* pars.Nstep(input(timerange',:)));
    spike_idx = find( sum( cif_chunk > repmat(isi,length(timerange),1),2),1); % index of 1st spike
    
    if isempty(spike_idx) % nothing happens in this chunk need to move to the next chunk!
        t = t+chunk_size;
        continue;
    end
    
    % Something spiked in this chunk - need to update history,coupling and
    % input terms
    n_idx = find(cif_chunk(spike_idx,:) > isi); % neurons that spiked
    spikes(t+spike_idx,n_idx) = (1>0);

    if (t+spike_idx == pars.maxt)
        break;
    end    
    
    for n = n_idx
        
        hist_range = (t+spike_idx+1):min(t+spike_idx+pars.Mhist,pars.maxt);
        % Update postspike term
        ps_update = pars.postspike_basis(1:length(hist_range),:)*pars.hweights.postspike(:,n);

        hist_term(hist_range,n) = hist_term(hist_range,n) + ps_update;
        input(hist_range,n) = input(hist_range,n) + ps_update;

        % Update coupling terms
    
        for n2=[1:n-1 n+1:pars.Nneurons]
            nrev = n;
            if (n > n2)
                nrev = nrev - 1;
            end
            cp_update = pars.coupling_basis(1:length(hist_range),:)*pars.hweights.coupling(:,(n2-1)*(pars.Nneurons-1)+nrev);  
            hist_term(hist_range,n2) = hist_term(hist_range,n2) + cp_update;
            input(hist_range,n2) = input(hist_range,n2) + cp_update;
        end
    end
    t = t+spike_idx;
end
    



