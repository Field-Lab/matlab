piece = '2011-01-19-0';

%% d03str
d03str_cones_path = '_snle_acquisition_2011-01-19-0_data003_data003-bayes-msf_15.00-2011-01-19-0_data003';


%% d05str
d05str_cones_path = '_snle_acquisition_2011-01-19-0_data005_data005-bayes-msf_15.00-2011-01-19-0_data005';


%% d09
d09_cones_path = '2011-01-19-0_data009_data009-bayes-msf_20.00-2011-01-19-0_data009_nwpca';


%% Stimuli

% C/R
stimuli = {'s06' 's07'};

% Wack
% Not clear which map was used :(
%stimuli(end+1) = {'s08'};

mapdirs.s06 = 'data003-off-midget-data006';
mapdirs.s07 = 'data005-off-midget-data007';
%mapdirs.s08 = 'data005-off-midget-data008';
for i = 1:length(stimuli)
    stim = read_stim_lisp_output([piece '/' stimuli{i}]);

    stim.mapdir = [server_data_path piece '/' mapdirs.(stimuli{i})];
    
    mapfiles = dir([stim.mapdir '/map-*.txt']);
    stim.mapfiles = {mapfiles.name};
    for j = 1:length(stim.mapfiles)
        stim.mapims{j} = dlmread([stim.mapdir '/' stim.mapfiles{j}]);
        fprintf('.');
    end
    fprintf('\n');
    
    assignin('caller', stimuli{i}, stim);
end; clear stim stimuli i j mapdirs mapfiles