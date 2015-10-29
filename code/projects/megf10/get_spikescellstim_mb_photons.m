function [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb_photons(datarun,cellids,timedur, bin_size)
%
% [NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb_photons(datarun,cellids,timedur, bin_size)
%
%
% Calculate average total number of spikes for each stimulus type      
%
% Input: Moving bar datarun structure
%
%   timedur: When last trial should stop
%        
%   cellids: cells being analyzed
%        
%   binsize: size within which to bin the data
%
% Outputs: NumSpikesCell: Has average total number of spikes for each stimulus type for each cell
%
%   StimComb: All the DG stimulus combinations (col 1: spatial,col 2: temporal, col 3: direction)
%
% Sneha Ravi 
% Last revision: 12-17-2012
% GDF: 2015-09-15         
          
%Initialize output matrices          
NumSpikesCell = zeros(length(cellids), length(datarun.stimulus.combinations));
MaxRate = zeros(length(cellids), length(datarun.stimulus.combinations));
StimComb = zeros(length(datarun.stimulus.combinations), length(fieldnames(datarun.stimulus.params)));
StimComb(:,1) = [datarun.stimulus.combinations.BAR_WIDTH];
StimComb(:,2) = [datarun.stimulus.combinations.DELTA];
StimComb(:,3) = [datarun.stimulus.combinations.DIRECTION];
%edges = 0:bin_size:500;

%Calculation of total spike number for each trial for each cell, then average calculated and placed in NumSpikesCell
temp_indices = get_cell_indices(datarun, cellids);
num_stimulus_combinations = length(datarun.stimulus.combinations);
reps = datarun.stimulus.repetitions;

date_check = input('is this data from before 2015-09-01?, y or n : ', 's');


if strcmp(date_check, 'y')
       % rgc loops over cells
    for rgc = 1:length(temp_indices)
        NumSpikesAll = zeros(1, length(datarun.stimulus.trials));
        MaxRateAll = zeros(1, length(datarun.stimulus.trials));
        % loop over epochs
        for epoch = 1:length(datarun.stimulus.triggers) 
            if (epoch == length(datarun.stimulus.triggers)) % handles the last trigger case
                % get spike indices between triggers
                Spikes_idx = datarun.spikes{rgc} >= datarun.stimulus.triggers(epoch) & datarun.spikes{rgc} < timedur;
                epoch_duration = timedur - datarun.stimulus.triggers(epoch);
                % get spike times
                Spikes_temp = datarun.spikes{rgc}(Spikes_idx);
                % subtract off time of first trigger so that start is zero 
                Spikes_temp = Spikes_temp - datarun.stimulus.triggers(epoch);
                % compute spike number
                NumSpikesAll(epoch) = sum(Spikes_idx);
                MaxRateAll(epoch) = max(histc(Spikes_temp, 0:bin_size:epoch_duration));
            else
                % get spike indices between triggers
                Spikes_idx = datarun.spikes{rgc} >= datarun.stimulus.triggers(epoch) & datarun.spikes{rgc} < datarun.stimulus.triggers(epoch+1);
                epoch_duration = datarun.stimulus.triggers(epoch+1) - datarun.stimulus.triggers(epoch);
                % get spike times
                Spikes_temp = datarun.spikes{rgc}(Spikes_idx);
                % subtract off time of first trigger so that start is zero 
                Spikes_temp = Spikes_temp - datarun.stimulus.triggers(epoch);
                % compute spike number
                NumSpikesAll(epoch) = sum(Spikes_idx);
                MaxRateAll(epoch) = max(histc(Spikes_temp, 0:bin_size:epoch_duration));
            end
        end
        size(NumSpikesAll)
        size(MaxRateAll)
        pause
        NumSpikesCell(rgc,:) = mean(reshape(NumSpikesAll, num_stimulus_combinations,reps),2);
        MaxRate(rgc,:) = mean(reshape(MaxRateAll, num_stimulus_combinations, reps),2);
    end 
    
    
    
    
    
else   
    % rgc loops over cells
    for rgc = 1:length(temp_indices)
        NumSpikesAll = zeros(1, length(datarun.stimulus.trials));
        MaxRateAll = zeros(1, length(datarun.stimulus.trials));
        trial = num_stimulus_combinations;
        pass = 0;
        % loop over epochs
        for epoch = 1:length(datarun.triggers) % the first trigger is garbage
            if trial == num_stimulus_combinations
                trial = 0;
                pass = pass + 1;
            elseif(epoch == length(datarun.triggers)) % handles the last trigger case
                % get spike indices between triggers
                Spikes_idx = datarun.spikes{rgc} >= datarun.triggers(epoch) & datarun.spikes{rgc} < timedur;
                epoch_duration = timedur - datarun.triggers(epoch);
                % get spike times
                Spikes_temp = datarun.spikes{rgc}(Spikes_idx);
                % subtract off time of first trigger so that start is zero 
                Spikes_temp = Spikes_temp - datarun.triggers(epoch);
                % compute spike number
                NumSpikesAll(trial + ((pass-1) * num_stimulus_combinations)) = sum(Spikes_idx);
                MaxRateAll(trial + ((pass-1) * num_stimulus_combinations)) = max(histc(Spikes_temp, 0:bin_size:epoch_duration));
            else
                trial = trial + 1;
                % get spike indices between triggers
                Spikes_idx = datarun.spikes{rgc} >= datarun.triggers(epoch) & datarun.spikes{rgc} < datarun.triggers(epoch+1);
                epoch_duration = datarun.triggers(epoch+1) - datarun.triggers(epoch);
                % get spike times
                Spikes_temp = datarun.spikes{rgc}(Spikes_idx);
                % subtract off time of first trigger so that start is zero 
                Spikes_temp = Spikes_temp - datarun.triggers(epoch);
                % compute spike number
                NumSpikesAll(trial + ((pass-1) * num_stimulus_combinations)) = sum(Spikes_idx);
                MaxRateAll(trial + ((pass-1) * num_stimulus_combinations)) = max(histc(Spikes_temp, 0:bin_size:epoch_duration));
            end
        end
        NumSpikesCell(rgc,:) = mean(reshape(NumSpikesAll, num_stimulus_combinations,reps),2);
        MaxRate(rgc,:) = mean(reshape(MaxRateAll, num_stimulus_combinations, reps),2);
    end
end


