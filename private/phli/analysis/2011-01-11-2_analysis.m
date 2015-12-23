piece = '2011-01-11-2';

%% d06: Voronoi contrast sensitivity:  6 contrasts (plus and minus) x 4 map
% files and blank trials, block-interleaved. 1 trial / sec, 2500 s (off midget)
d06 = load_data('2011-01-11-2/data006');
d06 = load_params(d06);
d06 = load_neurons(d06);
d06 = load_ei(d06, []);


%% d08: Voronoi contrast sensitivity:  6 contrasts (plus and minus) x 3 map
% files and blank trials, block-interleaved. 1 trial / sec, 1900 s (off midget)
d08 = load_data('2011-01-11-2/data008');
d08 = load_params(d08);
d08 = load_neurons(d08);
d08 = load_ei(d08, []);


%% d09: Voronoi contrast sensitivity:  6 contrasts (plus and minus) x 4 map files and blank trials, block-interleaved. 1 trial / sec, 1900 s (off parasol)
d09 = load_data('2011-01-11-2/data009');
d09 = load_params(d09);
d09 = load_neurons(d09);

% EIs need to be rerun!
d09 = load_ei(d09, []);


%% d10: Binary RGB 1-8-0.48-11111 ndf 0.0 3600 s
d10 = load_data('2011-01-11-2/data010');
d10 = load_params(d10);
d10 = load_neurons(d10);

% EIs need to be rerun!
d10 = load_ei(d10, []);


%% Stimuli

% C/R
stimuli = {'s06' 's08' 's09'};

mapdirs.s06 = 'data003-off-midget';
mapdirs.s08 = 'data005-off-midget';
mapdirs.s09 = 'data005-off-parasol';
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