
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, txt] = xlsread('/Users/colleen/Documents/Large Cell Data.xlsx', 'All large cell pieces');

for j= 2:size(txt,1)
    % piece = txt(j,1:3);
    run_opts.date=strtrim(txt{j,1}); % one slash at the end
    temp1 = strtrim(txt{j,2});
   temp2 =  strtrim(txt{j,3});
    run_opts.concatname=[temp1, temp2]; % Name (or modified name) of run, no slashes
    
    % Sometimes the data has two versions of the concate name
    run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
    % run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
    
    
    run_opts.save_location_root = '/Volumes/Lab/Users/crhoades/Cell Properties/';
    % Where to save
    run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, ];
    
    % run_opts.false_stixels = 0.5; % increase if stixels are size 32
    run_opts.num_frames = 30;
    run_opts.frames_past_zero = 0; %on the vision interface
    
    % Load Data
    datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
    datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
    datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
    opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
    try
        datarun=load_data(datarun,opt);
        
        run_opts.refresh = datarun.stimulus.interval/120*1000;
        % cell_indices = get_cell_indices(datarun,run_opts.cell_specification);
        
        for iter = 1:size(datarun.cell_types,2)
            run_opts.cell_specification = datarun.cell_types{iter}.name;
            run_opts.save_folder = [run_opts.date, '/', run_opts.concatname, '/', run_opts.cell_specification,  '/'];
            
            if ~exist([run_opts.save_location_root, run_opts.save_folder])
                mkdir([run_opts.save_location_root, run_opts.save_folder]);
            end
            
            cell_ids = get_cell_ids(datarun, run_opts.cell_specification);
            rf = zeros(length(cell_ids),1);
            t_zc= zeros(length(cell_ids),1);
            t_p = zeros(length(cell_ids),1);
            t_t = zeros(length(cell_ids),1);
            bi_ind = zeros(length(cell_ids),1);
            fr = zeros(length(cell_ids),1);
            amp = zeros(length(cell_ids),1);
            for i = 1:length(cell_ids)
                [rf(i), t_zc(i), t_p(i), t_t(i), bi_ind(i), fr(i), amp(i)] = get_timecourse_prop(datarun, cell_ids(i), run_opts);
            end
            
                     output.date = run_opts.date;
            output.concatname = run_opts.concatname;
            % if sum(strcmp(run_opts.cell_specification, {'ON parasol', 'ON midget', 'OFF midget', 'OFF parasol', 'OFF large 1','OFF large 2', 'OFF large 3', 'OFF large 4' ,'ON large 1', 'ON large 2', 'ON large 3', 'ON large 4'})) ~=1;
            %     if (strfind('ON', run_opts.cell_specification{1}) | strfind('on', run_opts.cell_specification{1}) | strfind('On', run_opts.cell_specification{1})) & (strfind('parasol', run_opts.cell_specification{1}) | strfind('Parasol', run_opts.cell_specification{1}))
            % end
            
            %     else
            
            output.cell_type = run_opts.cell_specification;
            output.parameters.rf = rf;
            output.parameters.t_zc = t_zc;
            output.parameters.t_p = t_p;
            output.parameters.t_t = t_t;
            output.parameters.bi_ind = bi_ind;
            output.parameters.fr = fr;
            
                        output.parameters.amp = amp;
            % end
            save([run_opts.save_location_root, run_opts.save_folder, 'output.mat'], 'output');
        end
    catch
        disp([run_opts.date, run_opts.concatname]);
    end
    clear datarun
end
