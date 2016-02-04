clear all
close all
% Load the Data for the DS run
run_opt.data_set = '2016-01-05-5';
run_opt.data_run = 'data003';
location = 'Data';

interval = 5; % number of seconds each stimulus was displayed for
repeats = 8;
% Where to save the data
location= ['/Volumes/Lab/Users/crhoades/DS_baden/', run_opt.data_set, '/', run_opt.data_run, '/'];

% You can give the cells as all of them (datarun.cell_ids) or give
% specific vision ids
% Find the cell to run by mapping a large cell EI from a white noise run

num_rand = 500;
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.params'];

% datarun.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', '08'];%run_opt.data_run(end-1:end)];
datarun.names.stimulus_path=['~/Desktop/2016-01-05-5/','s03','.txt'];


% datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.params'];
% datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.sta'];
%
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.neurons'];
% datarun.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', run_opt.data_run(end-1:end)];
opt=struct('verbose',1,'load_sta', 1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);

% datarun{1}=load_data(datarun{1},opt);
datarun=load_data(datarun,opt);

% cells = [124];%'all'; % 'all' or give a vector of vision cell ids

cells = get_cell_ids(datarun, 'all');


% datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));

% If the trigger interval in LabView is set long enought (~6 seconds for 5
% second stimuli), then this trigger_iti_thr should be fine.
datarun=load_stim(datarun,'correction_incomplet_run', 0, 'trigger_iti_thr', 0.01); % manually set threshold until got right number of triggers

if strcmp(cells, 'all')
    cells = datarun.cell_ids;
end
try
    spatial_tested = fliplr(datarun.stimulus.params.SPATIAL_PERIOD);
    temporal_tested = fliplr(datarun.stimulus.params.TEMPORAL_PERIOD);
catch
    datarun.stimulus.params.DIRECTION = datarun.stimulus.params.direction;
    datarun.stimulus.params.SPATIAL_PERIOD = datarun.stimulus.params.spatial_period;
    datarun.stimulus.params.TEMPORAL_PERIOD = datarun.stimulus.params.temporal_period;
    for i = 1:length(datarun.stimulus.trials)
        datarun.stimulus.trials(i).DIRECTION = datarun.stimulus.trials(i).direction;
        datarun.stimulus.trials(i).SPATIAL_PERIOD = datarun.stimulus.trials(i).spatial_period;
        datarun.stimulus.trials(i).TEMPORAL_PERIOD = datarun.stimulus.trials(i).temporal_period;
        
    end
    spatial_tested = fliplr(datarun.stimulus.params.SPATIAL_PERIOD);
    temporal_tested = fliplr(datarun.stimulus.params.TEMPORAL_PERIOD);
    direction_tested = fliplr(datarun.stimulus.params.DIRECTION);
    
end


tr=datarun.triggers; % all start triggers
trigger_marks=[];
counter = 1;
while counter< length(datarun.triggers)-1
    trigger_marks = [trigger_marks datarun.triggers(counter)];
    while datarun.triggers(counter+1) - datarun.triggers(counter) < 2;
        counter = counter+1;
        if counter> length(datarun.triggers)-1
            break;
        end
    end
    counter = counter+1;
end


tr = trigger_marks;

outlier_ds = cell(length(cells), 1);
outlier_os = cell(length(cells), 1);

for x = 1:length(cells)
    disp([num2str(x), ' out of ', num2str(length(cells))]);
    %     fprintf('.');
    cell_to_run = cells(x);
    %     cell_indices1=get_cell_indices(datarun{1}, [cell_to_run]);
    cell_indices2=get_cell_indices(datarun,[cell_to_run]);
    
    % Get spikes times for the specified cell
    n = datarun.spikes{cell_indices2(1)}';
    psth_r = [];
    for z = 1:size(tr,2)-1 % skip the first 12 trials due to contrast adapation
        ind = find(tr(z) == datarun.triggers);
        spatial = datarun.stimulus.trials(z).SPATIAL_PERIOD;
        temporal = datarun.stimulus.trials(z).TEMPORAL_PERIOD;
        direction = datarun.stimulus.trials(z).DIRECTION;
        %         rgb = datarun.stimulus.trials(z).RGB;
        
        num_triggers = interval/(temporal/120); % Number of times a band passed over one spot
        sub_spacing= datarun.triggers(ind:ind+num_triggers) - datarun.triggers(ind); % How to divide by the trial into the timing of the bands
        start = 0;
        stop_time = temporal/120;
        
        
        
        % first trial at the top
        for i=1:length(sub_spacing)-1,
            
            h=n-tr(z)-sub_spacing(i); % align the spike times
            stop = sub_spacing(i+1) - sub_spacing(i);
            if i == length(sub_spacing)-1
                stop = stop_time;
            end
            
            hh=find(h>=start & h<=stop); % Find spikes for each grating that passed over the spot
            psth_r=[psth_r; (h(hh))', repmat(length(sub_spacing)-i,[length(hh),1]);];
            %             psth=[psth_r; (h(hh))', repmat(length(sub_spacing)-i,[length(hh),1]);];
            
            
        end
        %         rasterplot(psth_r(:,1)+ (psth_r(:,2)-1)*stop_time,interval/stop_time ,stop_time)
        psth_mat{z} = psth_r; % Put all the psths from all the temporal, spatial and direction combinations into a cell array
        %                if rgb(1) == 0.48
        %                    psth_mat{z} = [];
        %                end
        
        psth_r = [];
        
    end
    
    
    %     for i = 25:size(tr,2)
    %             figure
    %             plot(psth_mat{i}(:,1), psth_mat{i}(:,2), 'o')
    % rasterplot(psth_r(:,1)+ (psth_r(:,2)-1)*0.25,20 ,0.25)
    %         end
    
    
    % Organize the data into [direction, spatial, temporal, number of spikes in
    % 5 second trial]
    by_trial = zeros(size(tr,2),4);
    for i = 1:size(tr,2)-1
        by_trial(i,1) = datarun.stimulus.trials(i).DIRECTION;
        by_trial(i,2) = datarun.stimulus.trials(i).SPATIAL_PERIOD;
        by_trial(i,3) = datarun.stimulus.trials(i).TEMPORAL_PERIOD;
        by_trial(i,4) = size(psth_mat{i},1); % number of spikes in that 5 seconds regardless of timing
        
    end
    
    by_trial_cell = cell(size(tr,2),4);
    for i = 1:size(tr,2)-1
        by_trial_cell{i,1} = datarun.stimulus.trials(i).DIRECTION;
        by_trial_cell{i,2} = datarun.stimulus.trials(i).SPATIAL_PERIOD;
        by_trial_cell{i,3} = datarun.stimulus.trials(i).TEMPORAL_PERIOD;
        by_trial_cell{i,4} = psth_mat{i}(:,1); % number of spikes in that 5 seconds regardless of timing
        
    end
    
    fig = figure;
    set(fig, 'PaperSize', [12 12], 'PaperPosition', [0 0 12, 12]);
    set(fig, 'Position', [0 0 1200 1200]);
    suptitle({[run_opt.data_set, ' ', run_opt.data_run];'DS'; ['Cell ', num2str(cells(x))]})
    
    
    fig_os = figure;
    set(fig_os, 'PaperSize', [12 12], 'PaperPosition', [0 0 12, 12]);
    set(fig_os, 'Position', [0 0 1200 1200]);
    suptitle({[run_opt.data_set, ' ', run_opt.data_run];'OS';['Cell ', num2str(cells(x))]})
    
    
    counter =1;
    for t = 1:length(temporal_tested)
        spat64 = find(by_trial(:,3) == temporal_tested(t));
        
        for s = 1:length(spatial_tested)
            spike_interest = cell(1,length(direction_tested));
            spike_temp = cell(1,length(direction_tested));
            spike_rand = cell(num_rand,1);
            for n = 1:num_rand
                spike_rand{n} = spike_temp;
            end
            
            spac = find(by_trial(:,2) == spatial_tested(s));
            
            
            for d= 1:length(direction_tested);
                dir = find(by_trial(:,1) == direction_tested(d));
                A = intersect(spat64, dir);
                
                C = intersect(A, spac);
                C= C(1:repeats-1);
                for c = 1:length(C)
                    spike_interest{:,d} = [spike_interest{:,d}; by_trial_cell{C(c),4}];
                end
                
            end
            
            %% null distribution
            %                 dir = find(by_trial(:,1) == direction_tested(d));
            %                 A = intersect(spat64, dir);
            
            C = intersect(spat64, spac);
            for n = 1:num_rand
                C_rand = C(randperm(length(C)));
                rand_C = reshape(C_rand(1:floor(size(C_rand(:),1)/length(direction_tested))*length(direction_tested)), length(direction_tested), floor(size(C_rand(:),1)/length(direction_tested)));
                for c = 1:size(rand_C,1)
                    spike_rand{n}{:,c} = [spike_rand{n}{:,c}; by_trial_cell{C_rand(c,:),4}];
                end
            end
            
            
            
            
            
            for d = 1:1:length(direction_tested);
                M(:,d) = hist(spike_interest{1,d});
            end
            
            % debugging
            %
            %     M =  20*ones(10,8);
            %     M(:,1) = 360;
            %      M(:,2) = 200;
            %      M(:,8) = 200;
            %           M(:,5) = 200;
            
            %     M = [ 270   187     6   105    96     4     1   274
            %    248    64     6    82   115    12    27   276
            %    322    51     6   199   215     5    26   304
            %    370     5    10    12     2   183   204   170
            %    288     3     1   172   312     6     0   240
            %    119     1     1    35   331    92    25   294
            %    259     4     3   277   280     8     3   213
            %    389     4     2   285   261     3     2    82
            %    309     4     8   380   140     4     1    10
            %    223    19     3     2     6   117   203   167];
            
            
            M_rand = cell(num_rand,1);
            for n = 1:num_rand
                for d = 1:1:length(direction_tested);
                    M_rand{n}(:,d) = hist(spike_rand{n}{1,d});
                end
            end
            
            
            [U,S,V] = svd(M);
            M2 = S(1,1)*U(:,1)*V(:,1)';
            
            clear i
            phi = exp(i*direction_tested*pi/180);
            phi_os = exp(2*i*direction_tested*pi/180);
            
            K = phi*V(:,1);
            K_os = phi_os*V(:,1);
            DSI{x}{s,t} = abs(K); % Yay!!
            OSI{x}{s,t} = abs(K_os); % Yay!!
            
            DSI_rand{x}{s,t} = nan(num_rand,1);
            OSI_rand{x}{s,t} = nan(num_rand,1);
            
            for n= 1:num_rand
                [U_rand,S_rand,V_rand] = svd(M_rand{n});
                M2_rand = S_rand(1,1)*U_rand(:,1)*V_rand(:,1)';
                
                clear i
                phi_rand = exp(i*direction_tested*pi/180);
                phi_rand_os = exp(2*i*direction_tested*pi/180);
                
                K_rand = phi_rand*V_rand(:,1);
                K_rand_os = phi_rand_os*V_rand(:,1);
                
                DSI_rand{x}{s,t}(n) = abs(K_rand); % Yay!!
                OSI_rand{x}{s,t}(n) = abs(K_rand_os); % Yay!!
                
            end
            
            DS_perc{x}{s,t} = prctile(DSI_rand{x}{s,t}, 95);
            OS_perc{x}{s,t} = prctile(OSI_rand{x}{s,t}, 95);
            
            if DS_perc{x}{s,t}<DSI{x}{s,t}
                outlier_ds{x}(s,t) = 1;
            else
                outlier_ds{x}(s,t) = 0;
            end
            
            if OS_perc{x}{s,t}<OSI{x}{s,t}
                outlier_os{x}(s,t) = 1;
            else
                outlier_os{x}(s,t) = 0;
            end
            
            figure(1)
            subplot(length(temporal_tested), length(spatial_tested), counter)
            
            
            hist(DSI_rand{x}{s,t});
            %             temp=get(fig,'Children')
            %             area = get(temp(1), 'Children');
            %             set(area, 'alpha', 0.2)
            %             set(child,'FaceAlpha',0.2)
            hold on
            y_axis = get(gca, 'ylim');
            x_axis = get(gca, 'xlim');
            plot([DSI{x}{s,t} DSI{x}{s,t}], y_axis, 'g', 'linewidth', 2)
            
            
            harea = patch([DS_perc{x}{s,t} DS_perc{x}{s,t} x_axis(2) x_axis(2)], [0 y_axis(2) y_axis(2) 0], 'r');
            set(harea, 'facealpha', 0.2)
            set(harea, 'Linestyle', 'none')
            hold on
            
            title(['Spatial: ',num2str(spatial_tested(s)) , ' / Temporal: ',num2str(temporal_tested(t))])
            
            xlabel('K')
            ylabel('Number of trials')
            
            
            figure(2)
            subplot(length(temporal_tested), length(spatial_tested), counter)
            hist(OSI_rand{x}{s,t});
            hold on
            y_axis = get(gca, 'ylim');
            x_axis = get(gca, 'xlim');
            
            plot([OSI{x}{s,t} OSI{x}{s,t}], y_axis, 'g')
            
            harea = patch([OS_perc{x}{s,t} OS_perc{x}{s,t} x_axis(2) x_axis(2)], [0 y_axis(2) y_axis(2) 0], 'r');
            set(harea, 'facealpha', 0.2)
            set(harea, 'Linestyle', 'none')
            hold on
            
            title(['Spatial: ',num2str(spatial_tested(s)) , ' / Temporal: ',num2str(temporal_tested(t))])
            
            xlabel('K')
            ylabel('Number of trials')
            
            counter = counter+1;
        end
    end
    
    filepath_DS = [location, '/DS/Cell_', num2str(cells(x))];
    filepath_OS = [location, '/OS/Cell_', num2str(cells(x))];
    
    if ~exist([location, '/DS'])
        mkdir([location, '/DS'])
    end
    
    if ~exist([location, '/OS'])
        mkdir([location, '/OS'])
    end
    
    % set(gcf, 'renderer', 'opengl')
    
    %     export_fig([filepath, 'Cell_',num2str(cell_to_run)], '-pdf')
    figure(1)
    print(fig,'-dpdf',[filepath_DS]);
    figure(2)
    print(fig_os,'-dpdf',[filepath_OS]);
    
    
    
    close all
end

save([location, '/DS/all'], 'outlier_ds');
save([location, '/OS/all'], 'outlier_os');


%%

outlier_ds_mat = cell2mat(outlier_ds);

answer(:,1:3) = outlier_ds_mat(1:2:end, :);
answer(:,4:6) = outlier_ds_mat(2:2:end, :);

 [Y,I] = sortrows(answer);

 cells = get_cell_ids(datarun, 'all');

  cells(I(size(I,1)-12:size(I,1)))
  
  
  outlier_os_mat = cell2mat(outlier_os);

answer(:,1:3) = outlier_os_mat(1:2:end, :);
answer(:,4:6) = outlier_os_mat(2:2:end, :);

 [Y,I] = sortrows(answer);

 cells = get_cell_ids(datarun, 'all');

  cells(I(size(I,1)-13:size(I,1)))
  
