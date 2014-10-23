% DATA PARAMETERS
run_opt.load = true; % T/F
run_opt.data_set = '2007-08-24-4'; % '2007-03-27-1'
run_opt.data_run = 7; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
run_opt.config_num = 1; % 1-4
run_opt.cell_type = 'On midget'; % on/off parasol, on/off midget
run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
run_opt.auto_set = true; % T/F -- note: overwrites run_opt params

% NUMERICAL PARAMETERS
run_opt.tau = .01; % tuning parameter
run_opt.tol = 1e-3;
run_opt.trial_estimate_start = 120;
run_opt.velocity_lim = 150; % >0

% ANALYSES TO RUN
run_opt.downsample_spikes = false;
run_opt.raster = false; % T/F
run_opt.trial_estimate = false; % T/F


tic;

% Auto set parameters
if run_opt.auto_set
    [run_opt.cell_types, run_opt.velocity_lim, run_opt.config_num, run_opt.trial_estimate_start, run_opt.tol] =...
        auto_set_params(run_opt.data_set, run_opt.data_run);
end

% Load data fresh
if run_opt.load

    clear datarun tr

    if strcmp(run_opt.data_set, '2007-03-27-1')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-03-27-1/data011-nwpca/data011-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-08-24-4/data001-nwpca/data001-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-08-24-4/data%03d-from-data001-nwpca/data%03d-from-data001-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-08-24-4/Stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2005-04-26-0')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2005-04-26-0/data';
    end
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
    datarun=load_data(datarun,opt);
    datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true)); 
    datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0); 
    
end

if run_opt.raster || run_opt.trial_raster || run_opt.trial_estimate
    
    % Get indices for specified cell type and order by RF position
    cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
    cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
    
    cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits);
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1));
    
    cell_indices1 = cell_indices1(cell_sort_idx);
    cell_indices2 = cell_indices2(cell_sort_idx);
    
    % Find trial start and stop times
    start = 0;
    stop = mean(datarun{2}.triggers(2:2:end) - datarun{2}.triggers(1:2:end));

    tr=datarun{2}.triggers(1:2:end); % triggers mark the beginning and end
    t=find(datarun{2}.stimulus.trial_list==run_opt.config_num);
    tr=tr(t);
    
end

% downsample spikes
if run_opt.downsample_spikes
    if strcmp(run_opt.data_set, '2007-03-27-1')
        n1 = round(mean_spikes_onp_exp1(run_opt.data_run-13));
        n2 = round(mean_spikes_onm_exp1(run_opt.data_run-13));
        if n1>n2
            if strcmp(run_opt.cell_type, 'On parasol')
                for i=cell_indices2
                    % subtract difference in means:
                    %datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},length(datarun{2}.spikes{i})-abs(n1-n2),'Replace',false);
                    % fractional downsampling:
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},length(datarun{2}.spikes{i})*n2/n1,'Replace',false);
                end
            else
                run_opt.trial_estimate = false;
            end
        elseif n2>n1
            if strcmp(run_opt.cell_type, 'On midget')
                for i=cell_indices2
                    % subtract difference in means:
                    %datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},length(datarun{2}.spikes{i})-abs(n1-n2),'Replace',false);
                    % fractional downsampling
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},length(datarun{2}.spikes{i})*n1/n2,'Replace',false);                    
                end
            else
                run_opt.trial_estimate = false;
            end
        end
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        n1 = round(mean_spikes_onp_exp2(run_opt.data_run-3));
        n2 = round(mean_spikes_onm_exp2(run_opt.data_run-3));
        if n1>n2
            if strcmp(run_opt.cell_type, 'On parasol')
                for i=cell_indices2
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},max(length(datarun{2}.spikes{i})-abs(n1-n2),0),'Replace',false);
                end
            else
                run_opt.trial_estimate = false;
            end
        elseif n2>n1    
            if strcmp(run_opt.cell_type, 'On midget')
                for i=cell_indices2
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},length(datarun{2}.spikes{i})-abs(n1-n2),'Replace',false);
                end
            else
                run_opt.trial_estimate = false;
            end
        end
    end
end

if run_opt.raster %raster

    k=1; kmin=1; kmax=length(cell_indices2); hk=loop_slider(k,kmin,kmax);

    while k
        if ~ishandle(hk)
            break
        end
        k=round(get(hk,'Value')); 
 
        psth_raster(start,stop,datarun{2}.spikes{cell_indices2(k)}',tr);
        title(sprintf('%d %.2f', datarun{2}.cell_ids(cell_indices2(k)), datarun{1}.vision.sta_fits{cell_indices1(k)}.mean(1) ))

        uiwait;
    end
end

if run_opt.trial_estimate
    
    % start parallel pool
    poolobj = parpool;
    
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 30, 'LargeScale', 'off');
    estimates = zeros(size(tr));
    spikes = datarun{2}.spikes;
    parfor i = 1:length(tr)
        estimates(i) = fminunc(@(v) -pop_motion_signal(v, spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
        fprintf('for trial %d, the estimated speed was %d', i, estimates(i))
    end
    
    % stop parallel pool
    delete(poolobj);
    
    % save estimates
    save(sprintf('/home/vision/Dropbox/Lab/Development/matl ab-standard/private/malcolm/results/%s/downsample_%s_data_run_%02d_config_%d.mat', run_opt.data_set, run_opt.cell_type, run_opt.data_run, run_opt.config_num), 'estimates')
end

ElapsedTime=toc;