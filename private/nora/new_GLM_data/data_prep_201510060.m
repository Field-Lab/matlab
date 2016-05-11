clear
Analysis_Path = '/Volumes/Analysis/2015-10-06-0/data000-data015-norefit/';
class = 'data000';
NSEM_runs = {'data002', 'data007', 'data008', 'data009', 'data012', 'data015'};
WN_runs = {'data003', 'data010', 'data013'};
WN_seeds = [11111 22222 33333];
WN_fit_frames = [3600 6000 3600];
extra_string = '-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015';
datarun_class = load_data([Analysis_Path class extra_string '/' class extra_string], struct('load_neurons', 0, 'load_params', 1));
reps = 60; 
idx = 1:reps;
for cell_type = {'On Parasol'}
    cells = get_cell_ids(datarun_class, cell_type{1}); % cell ids to fit
    datarun_class = load_sta(datarun_class, 'load_sta', cells);
    for run = 1:3
        datarun = load_data([Analysis_Path WN_runs{run} extra_string '/' WN_runs{run} extra_string], struct('load_neurons', 1, 'load_params', 1));
        prepped_data = interleaved_data_prep(datarun, [WN_fit_frames(run) 1200], reps,'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1','seed', 11111, ...
            'variable_seed_start', 11111*run, 'STA_check', 0);
        fitspikes_cell(idx, :) = prepped_data.fitspikes;
        fitmovie_cell(idx) = prepped_data.fitmovie;
        testmovie = prepped_data.testmovie;
        monitor_refresh(run) = 100/median(diff(datarun.triggers));
        idx = idx+reps;
    end
end

%
