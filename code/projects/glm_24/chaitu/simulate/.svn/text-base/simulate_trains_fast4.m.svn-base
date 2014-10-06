function [exflag spikes input hist_term] = simulate_trains_fast4(pars,dt,kx,b,Deff,ps_coeffs,cp_coeffs,cpl_mat,chunk_size)

% Function for simulating spike trains given a filtered stimulus
% Chaitu Ekanadham 07/27/2009

% Args: pars - model parameters
%       dt - timestep
%       kx - filtered stimulus (time x neuron)
%       Deff - spikes of neurons not involved in the simulation - these are considered given data like the stimulus
%       ps_coeffs - postspike basis coefficients (pars.nofilters_postspike x neuron)
%       cp_coeffs - cell array with coupling basis coefficients (nth element is pars.nofilters_coupling x nc-1 where nc is the number of 1s  in the nth row of cpl_mat)
%       chunk size - for time-rescaling procedure
%       cpl_mat - logical matrix specifying which neuron pairs have a coupling filter

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
DEFAULT_CHUNK_SIZE = 100;%max(100,floor(spT/10));
if (~exist('chunk_size','var') || chunk_size < 0)
    chunk_size = DEFAULT_CHUNK_SIZE;
end

%chunk_size
t = 0; % time index (not ACTUAL time!)

% Run simulation till end of time interval

exflag = 1;

Neff = size(spikes,2);

PS = pars.postspike_basis * ps_coeffs;
%if (Neff > 1)
%    CP = pars.coupling_basis * cp_coeffs;
%end

counter = 0;
1;
fprintf('\n');
disppercent(-inf,'Starting trial.');

while (t < spT)

    disppercent(t/spT);
    
    timerange = (t+1):min(t+chunk_size,spT);
    
    counter = counter + 1;
    
    %if (mod(counter,100) == 0)
    %    fprintf('Iter %d: spT-t=%f\n',counter,spT-t);
    %end
    
    if (max(max(input(timerange',:))) > 40)
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
    
    if (isempty(fixedspike_idx) || (~isempty(spike_idx) && spike_idx < fixedspike_idx)) % a real neuron spiked
        n_idx = find(cif_chunk(spike_idx,:) > isi); % index of neuron(s) that spiked
        spikes(t+spike_idx,n_idx) = true;
    else % a virtual (conditioned) neuron spiked
        n_idx = find(fixedspike_chunk(fixedspike_idx,:));
        spike_idx = fixedspike_idx;
    end
    
    

    %if (t+spike_idx == spT)
    %    fprintf('Exiting out of loop.\n');
    %    break;
    %end    
    
    for r = 1:length(n_idx) % Iterate through spiking neurons
        n = n_idx(r); % index of a neuron that spiked 
        hist_range = (t+spike_idx+1):min(t+spike_idx+pars.Mhist,spT); % time range where history update will take place
        
        if (n <= pars.Nneurons) % check if this is a real neuron (not a virtual)
                % Update postspike term
                hist_term(hist_range,n) = hist_term(hist_range,n) + PS(1:length(hist_range),n);
                input(hist_range,n) = input(hist_range,n) + PS(1:length(hist_range),n);
        end

        % Update coupling terms of simulated neurons
        if (Neff > 1)
            n2_idx = find(cpl_mat(:,n)); % neurons potentially effected by neuron n's spike by coupling connections
            
            for l=1:length(n2_idx)% iterate through potentially affected neurons [1:min(n-1,pars.Nneurons) n+1:pars.Nneurons]
                n2 = n2_idx(l);
                %nrev = n;
                %if (n > n2)
                %    nrev = nrev - 1;
                %end
                % Find the index of the coupling filter coeffs corresponding to
                % n's effect on n2
                n2neighbors = find(cpl_mat(n2,:));
                1;
                n2_n_coeffs = cp_coeffs{n2}(:,find(n2neighbors == n,1));
                cp_update = pars.coupling_basis(1:length(hist_range),:)*n2_n_coeffs;      %(:,(n2-1)*(Neff-1)+nrev);
                hist_term(hist_range,n2) = hist_term(hist_range,n2) + cp_update;
                input(hist_range,n2) = input(hist_range,n2) + cp_update;
            end
        end
    end
    t = t+spike_idx;
end
elapsedTime = disppercent(inf);
1;
exflag = 0;