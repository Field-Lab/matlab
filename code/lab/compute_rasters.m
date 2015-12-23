function datarun = compute_rasters(datarun, cell_spec, varargin)

%sneha ravi 09/12/12 changed datarun.triggers to datarun.stimulus.triggers
%in 2 instances

% Initialize the inputParser
p = inputParser;

% variable opts for wrapper function
p.addParamValue('triggers_per_trial', 1, @isnumeric);
p.addParamValue('compute_psth', true, @islogical);
p.addParamValue('print_summaries', false, @islogical);
p.addParamValue('plot_psth', true, @islogical);
p.addParamValue('foa', 1)
p.addParamValue('master_cell_ids', []);
p.addParamValue('save_path', '/Analysis/raster-files/');


% opts for rasters
p.addParamValue('start', 0, @isnumeric);
p.addParamValue('stop', [], @isnumeric);
p.addParamValue('tic_color', [0 0 0]);
p.addParamValue('line_space', 0.2, @isnumeric)
p.addParamValue('tic_thickness', 1);
p.addParamValue('first_trial', 1, @isnumeric)

% more opts for psths
p.addParamValue('bin_size', 0.1, @isnumeric);
p.addParamValue('hist_color', [0 0 0]);

% parse inputs
p.parse(varargin{:});
params = p.Results;

raster_opts.start = params.start;
raster_opts.stop = params.stop;
raster_opts.tic_color = params.tic_color;
raster_opts.line_space = params.line_space;
raster_opts.tic_thickness = params.tic_thickness;
raster_opts.first_trial = params.first_trial;
raster_opts.plot = 'true';
raster_opts.foa = -1;

psth_opts.start = params.start;
psth_opts.stop = params.stop;
psth_opts.first_trial = params.first_trial;
psth_opts.bin_size = params.bin_size;
psth_opts.hist_color = params.hist_color;
psth_opts.plot_hist = params.plot_psth;

psth_opts.foa = -1;


% check to see if trial structure already exists
if ~isfield(datarun, 'trials');
    datarun.trials = [];
end    

rasters = cell(length(datarun.cell_ids), 1);
psths = cell(length(datarun.cell_ids), 1);


%%%% begin function %%%%

cell_indices = get_cell_indices(datarun, cell_spec);

trig_indices = 1:params.triggers_per_trial:length(datarun.stimulus.triggers);
begin_trial_triggers = datarun.stimulus.triggers(trig_indices);

for rgc = 1:length(cell_indices);
    figure(1); clf; subplot(2,1,1)
    temp_spike_times = datarun.spikes{cell_indices(rgc)};
    rgc_raster = get_raster(temp_spike_times, begin_trial_triggers, raster_opts);
    rasters{cell_indices(rgc)} = rgc_raster;
    
    subplot(2,1,2)
    if params.compute_psth
        rgc_psth = get_psth(temp_spike_times, begin_trial_triggers, psth_opts);
    end
    psths{cell_indices(rgc)} = rgc_psth;
    
    subplot(2,1,1)
    

    if params.print_summaries
        if ~isempty(params.master_cell_ids);
            file_name = [num2str(params.master_cell_ids(rgc)),'.pdf'];
            title(['cell ID ', num2str(datarun.cell_ids(cell_indices(rgc))),' -> ',num2str(params.master_cell_ids(rgc))])
        else
            file_name = num2str(datarun.cell_ids(cell_indices(rgc)));
            title(['cell ID ', num2str(datarun.cell_ids(cell_indices(rgc)))])
        end
        print(1, [params.save_path,file_name], '-dpdf'); 
    end
    
end

if ~isfield(datarun.trials, 'rasters')
    datarun.trials.rasters = rasters;
else
    for cnt = 1:length(cell_indices)
        datarun.trials.rasters{cell_indices(cnt)} = rasters{cell_indices(cnt)};
    end
end

if ~isfield(datarun.trials, 'psths');
    datarun.trials.psths = psths;
else
    for cnt = 1:length(cell_indices)
        datarun.trials.psths{cell_indices(cnt)} = psths{cell_indices(cnt)};
    end
end


    
    
    
    
    
