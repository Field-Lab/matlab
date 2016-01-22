function [exflag spikes input hist_term] = simulate_trains_fast3(pars,dt,kx,b,Deff,ps_coeffs,cp_coeffs,chunk_size)

% Function for simulating spike trains given a filtered stimulus
% Args: model parameters
%       timestep dt
%       filtered stimulus x
%       Deff - spikes of neurons not involved in the simulation - these are considered given data like the stimulus
%       ps_coeffs - postspike basis coefficients (pars.nofilters_postspike x neuron)
%       cp_coeffs - coupling basis coefficients (pars.nofilters_coupling x (neuron1*neuron2))
%       chunk size for time-rescaling procedure
% Returns:
% spikes - binary matrix (time x neuron) of spikes 
% input - log of the CIF (time x neuron)
% hist_term - h currents (time x neuron) (for debugging mostly)

% Initialize the simulation variables
spT = pars.maxt*pars.fac;
spikes = sparse(logical(false(spT,pars.Nneurons))); % the spike trains (logical matrix)
spikes = [spikes Deff]; % attach the other neuron spikes
input = repmat(b',spT,1) + reprows(kx,pars.fac); % the input current (log of the CIF) pars.maxt x pars.Nneurons - can interpolate here if you want
hist_term = zeros(spT,pars.Nneurons); % the history terms
%DEBUGGING
%hist_test = zeros(pars.maxt,Ntrials,pars.Nneurons);




% Cap if necessary (only do this for the first time step - there is an
% error if we have to do this in the middle)
cap_idx = input(1,:) > log(1/dt);
input(1,cap_idx) = log(1/dt);
%if (sum(sum(sum(cap_idx))) > 0)
%    fprintf('Capped the initial CIF value for %d neurons to %d\n\n',sum(sum(sum(cap_idx))),1/dt);
%end

% Time window to stimulate per iteration
DEFAULT_CHUNK_SIZE = max(100,floor(spT/10));
if (~exist('chunk_size','var') || chunk_size < 0)
    chunk_size = DEFAULT_CHUNK_SIZE;
end

%chunk_size
t = 0; % time index (not ACTUAL time!)

% Run simulation till end of time interval

exflag = 1;

Neff = size(spikes,2);

PS = pars.postspike_basis * ps_coeffs;
if (Neff > 1)
    CP = pars.coupling_basis * cp_coeffs;
end

counter = 0;
1;

%disppercent(-inf,'Starting trial.');

while (t < spT)

%    disppercent(t/spT);
    
    timerange = (t+1):min(t+chunk_size,spT);
    
    counter = counter + 1;
    
    %if (mod(counter,100) == 0)
    %    fprintf('Iter %d: spT-t=%f\n',counter,spT-t);
    %end
    
    if (max(max(input(timerange',:))) > 20)
        1;
        return;
    end
    
    %t/spT
    
    % Generate a random Exp(1) value
    isi = -log(1-rand(1,pars.Nneurons));%exprnd(1,1,pars.Nneurons);
    % Find the first index at which the integrated cifs exceed the isi    
    cif_chunk = cumsum(dt .* pars.Nstep(input(timerange',:)));
    if (~isempty(Deff))
        fixedspike_chunk = Deff(timerange,:);
    else
        fixedspike_chunk = [];
    end
    
    1;
    spike_idx = find( sum( cif_chunk > repmat(isi,length(timerange),1),2),1); % index of 1st spike of simulated neurons
    fixedspike_idx = find(sum(fixedspike_chunk,2),1);
    
    if (isempty(spike_idx) && isempty(fixedspike_idx)) % nothing happens in this chunk need to move to the next chunk!
        t = t+chunk_size;
        continue;
    end
    
    if (isempty(spike_idx))
        1;
    end
    
    % Something spiked in this chunk - need to update history,coupling and
    % input terms
    1;
   
    %if (isempty(fixedspike_idx) || isempty(spike_idx))
    %    fixedspike_idx
    %    spike_idx
    %end
    
    if (isempty(fixedspike_idx) || (~isempty(spike_idx) && spike_idx < fixedspike_idx))
        n_idx = find(cif_chunk(spike_idx,:) > isi); % neurons that spiked
        spikes(t+spike_idx,n_idx) = true;
    else
        n_idx = find(fixedspike_chunk(fixedspike_idx,:));
        spike_idx = fixedspike_idx;
    end
    
    

    %if (t+spike_idx == spT)
    %    fprintf('Exiting out of loop.\n');
    %    break;
    %end    
    
    for n = n_idx
        hist_range = (t+spike_idx+1):min(t+spike_idx+pars.Mhist,spT);
        if (n <= pars.Nneurons)
                % Update postspike term
                hist_term(hist_range,n) = hist_term(hist_range,n) + PS(1:length(hist_range),n);
                input(hist_range,n) = input(hist_range,n) + PS(1:length(hist_range),n);
        end

        % Update coupling terms of simulated neurons
    
        for n2=[1:min(n-1,pars.Nneurons) n+1:pars.Nneurons]
            nrev = n;
            if (n > n2)
                nrev = nrev - 1;
            end
            cp_update = pars.coupling_basis(1:length(hist_range),:)*cp_coeffs(:,(n2-1)*(Neff-1)+nrev);  
            hist_term(hist_range,n2) = hist_term(hist_range,n2) + cp_update;
            input(hist_range,n2) = input(hist_range,n2) + cp_update;
        end
    end
    t = t+spike_idx;
end
%elapsedTime = disppercent(inf);
1;
exflag = 0;