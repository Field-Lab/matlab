
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[~, txt,array_size] = xlsread('/Volumes/Lab/Users/crhoades/Large Cell Data ARVO.xlsx');
array_size = cell2mat(array_size(2:end,4));
% txt = {'
% };

for j= 2:size(txt,1)


    % piece = txt(j,1:3);
    run_opts.date=strtrim(txt{j,1}); % one slash at the end
    temp1 = strtrim(txt{j,2});
    temp2 =  strtrim(txt{j,3});
    run_opts.concatname=temp1; % Name (or modified name) of run, no slashes\
    run_opts.dataname = temp2;
    
    % Sometimes the data has two versions of the concate name
    run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',run_opts.dataname, '/',run_opts.dataname];
    % run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
    
    
    run_opts.save_location_root = '/Volumes/Lab/Users/crhoades/Cell Properties/';
    % Where to save
    run_opts.filepath= [run_opts.save_location_root, run_opts.date, '/', run_opts.concatname, '/', run_opts.dataname];
    
    % run_opts.false_stixels = 0.5; % increase if stixels are size 32
    run_opts.num_frames = 30;
    run_opts.frames_past_zero = 0; %on the vision interface
    
    % Load Data
    datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opts.file_name, '.neurons'];
    datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];
    datarun.names.rrs_sta_path = ['/Volumes/Analysis/', run_opts.file_name, '.sta'];
    datarun.names.rrs_ei_path = ['/Volumes/Analysis/', run_opts.file_name, '.ei'];
    
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
    opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
    try
        datarun=load_data(datarun,opt);
        run_opts.refresh = datarun.stimulus.interval/120*1000;
        % cell_indices = get_cell_indices(datarun,run_opts.cell_specification);
        pca_data = [];
        for iter = 1:size(datarun.cell_types,2)
            clear probabilities;
            run_opts.cell_specification = datarun.cell_types{iter}.name;
            run_opts.save_folder = [run_opts.date, '/', run_opts.concatname, '/', run_opts.cell_specification,  '/'];
            
            if ~exist([run_opts.filepath,'/', run_opts.cell_specification]);
                
                mkdir([run_opts.filepath,'/', run_opts.cell_specification]);
            end
            
            cell_ids = get_cell_ids(datarun, run_opts.cell_specification);
            cell_indices = get_cell_indices(datarun, run_opts.cell_specification);
            
            
            rf = zeros(length(cell_ids),1);
            t_zc= zeros(length(cell_ids),1);
            t_p = zeros(length(cell_ids),1);
            t_t = zeros(length(cell_ids),1);
            bi_ind = zeros(length(cell_ids),1);
            fr = zeros(length(cell_ids),1);
            amp = zeros(length(cell_ids),1);
%             width = zeros(length(cell_ids),1);
%             isi_max = zeros(length(cell_ids),1);
            acf_mean=zeros(length(cell_ids),1);
            acf_var=zeros(length(cell_ids),1);
            acf_sd=zeros(length(cell_ids),1);
            acf_end_val=zeros(length(cell_ids),1);
            acf_peak=zeros(length(cell_ids),1);
            acf_peak_time=zeros(length(cell_ids),1);
            acf_rise=zeros(length(cell_ids),1);
            acf_slope_up = zeros(length(cell_ids),1);
            acf_slope_down = zeros(length(cell_ids),1);
            maximum= zeros(length(cell_ids),1);
            minimum= zeros(length(cell_ids),1);
            minimum_ind= zeros(length(cell_ids),1);
            maximum_ind= zeros(length(cell_ids),1);
            zc= zeros(length(cell_ids),1);
            bi_ind_ei= zeros(length(cell_ids),1);
            lobe1= zeros(length(cell_ids),1);
            lobe2= zeros(length(cell_ids),1);
            ratio_to_neighbors= zeros(length(cell_ids),1);
            
            norm= cell(length(cell_ids),1);
            cell_counter  = 0;
            for i = 1:length(cell_ids)
                
                try
                    [rf(i), t_zc(i), t_p(i), t_t(i), bi_ind(i), fr(i), amp(i)] = get_timecourse_prop(datarun, cell_ids(i), run_opts);
%                     [maximum(i), minimum(i), minimum_ind(i),  maximum_ind(i), zc(i), bi_ind_ei(i), lobe1(i), lobe2(i), ratio_to_neighbors(i)] = ei_properties(datarun, cell_ids(i), array_size(j));
                    
                    
                catch
                    disp([num2str(iter) '/' num2str(i)])
                end

                %                 try
                
                %                     [probabilities{i}, bins{i}, norm{i}] = autocorrelation(datarun.spikes{cell_indices(i)},0.001, 0.1, datarun.duration);
                %                     [width(i)] = inter_spike_interval(datarun.spikes{cell_indices(i)});
                %                     cell_counter = cell_counter +1;
                %
                %                     acf_mean(i) = norm{i}.mean;
                %                     acf_var(i) = norm{i}.var;
                %                     acf_sd(i) = norm{i}.sd;
                %                     acf_end_val(i) = norm{i}.end_value;
                %                     acf_peak(i) = norm{i}.peak;
                %                     acf_peak_time(i) = norm{i}.peak_time;
                %                     acf_rise(i) = norm{i}.rise;
                %                     acf_slope_up(i) = norm{i}.slope_up;
                %                     acf_slope_down(i) = norm{i}.slope_down;
                %                 catch
                probabilities{i} = nan(1,101);
                bins{i} = nan(1,101);
                cell_counter = cell_counter +1;
                width(i) =  nan;
                acf_mean(i) = nan;
                acf_var(i) = nan;
                acf_sd(i) = nan;
                acf_end_val(i) = nan;
                acf_peak(i) = nan;
                acf_peak_time(i) = nan;
                acf_rise(i) = nan;
                acf_slope_up(i) =nan;
                acf_slope_down(i) =nan;
                
                
                %                 end
                
            end
%             [~, elimin1] = removeoutliers(minimum_ind);
%             [~, elimin2] = removeoutliers(maximum_ind);
%             [~, elimin3] = removeoutliers(zc);
%             [~, elimin4] = removeoutliers(bi_ind);
%             remove = unique([elimin1; elimin2; elimin3; elimin4]);
%             
%             rf(remove) = nan;
%             t_zc(remove) = nan;
%             t_p(remove) = nan;
%            t_t(remove) = nan;
%             bi_ind(remove) = nan;
%             fr(remove) = nan;
%             amp(remove) = nan;
% %             output.parameters.width = width;
% %             output.parameters.isi_max = isi_max;
%             acf_mean(remove) = nan;
%             acf_var(remove) = nan;
%             acf_sd(remove) = nan;
%            acf_end_val(remove) = nan;
%             acf_peak(remove) = nan;
%            acf_peak_time(remove) = nan;
%             acf_rise(remove) = nan;
%             acf_slope_up(remove) = nan;
%             acf_slope_down(remove) = nan;
%             maximum(remove) = nan;
%             minimum(remove) = nan;
%             minimum_ind(remove) = nan;
%             maximum_ind(remove) = nan;
%             zc(remove) = nan; 
%             bi_ind_ei(remove) = nan;
%             lobe1(remove) = nan; 
%             lobe2(remove) = nan; 
%             ratio_to_neighbors(remove) = nan;
%             maximum(i), minimum(i), minimum_ind(i),  maximum_ind(i), zc(i), bi_ind(i), lobe1(i), lobe2(i), ratio_to_neighbors(i)
% maximum = removeoutliers(maximum);
% maximum = removeoutliers(maximum);
% maximum = removeoutliers(maximum);
%             count(iter) = cell_counter;
            %             try
            %                 test = cell2mat(probabilities);
            %                 test2 = reshape(test, length(bins{1}), length(test)/length(bins{1}));
            %
            %                 pca_data = [pca_data, test2];
            %             catch
            %                 disp('pca catch')
            %             end
            
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
%             output.parameters.width = width;
%             output.parameters.isi_max = isi_max;
            output.parameters.acf_mean =acf_mean;
            output.parameters.acf_var =acf_var;
            output.parameters.acf_sd =acf_sd;
            output.parameters.acf_end_val =acf_end_val;
            output.parameters.acf_peak =acf_peak;
            output.parameters.acf_peak_time =acf_peak_time;
            output.parameters.acf_rise =acf_rise;
            output.parameters.acf_slope_up =acf_slope_up;
            output.parameters.acf_slope_down =acf_slope_down;
            output.parameters.maximum = maximum;
            output.parameters.minimum =minimum;
            output.parameters.minimum_ind = minimum_ind;
            output.parameters.maximum_ind = maximum_ind;
            output.parameters.zc = zc;
            output.parameters.bi_ind_ei = bi_ind_ei;
            output.parameters.lobe1 = lobe1;
            output.parameters.lobe2 = lobe2;
            output.parameters.ratio_to_neighbors = ratio_to_neighbors;
            
            % end
            save([run_opts.filepath,'/', run_opts.cell_specification, '/', 'output.mat'], 'output');
        end
    catch
        disp([run_opts.date, run_opts.concatname]);
    end
    %     [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(pca_data', 'NumComponents',2);
    %     counter = 1;
    %     for iter = 1:size(datarun.cell_types,2)
    %         load([run_opts.filepath,'/',  datarun.cell_types{iter}.name, '/', 'output.mat']);
    %         cell_ids = get_cell_ids(datarun, datarun.cell_types{iter}.name);
    %
    %         output.parameters.pca1 = SCORE(counter:counter +count(iter) -1,1);
    %         output.parameters.pca2 = SCORE(counter:counter +count(iter) -1,2);
    %         save([run_opts.filepath,'/', datarun.cell_types{iter}.name, '/', 'output.mat'], 'output');
    %         counter = counter + count(iter);
    %     end
    
    clear datarun
end
