% runDS
% Colleen Rhoades
% June 2015
% rhoades@stanford.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of this file is to analyze the results from drifting gratings.

% You will see a graph showing the number of total spikes during the stimulus presentation versus the direction of the stimulus.
% If a direction is preferred by the cell, then the spike count should be higher.
% The plot is divided for different spatial and temporal frequencies.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the Data for the DS run
run_opt.data_set = '2016-01-05-5';
run_opt.data_run = 'data003';
location = 'Data';

interval = 5; % number of seconds each stimulus was displayed for

% Where to save the data
filepath= ['/Volumes/Lab/Users/crhoades/DS/', run_opt.data_set, '/', run_opt.data_run, '/'];

% You can give the cells as all of them (datarun.cell_ids) or give
% specific vision ids
% Find the cell to run by mapping a large cell EI from a white noise run
cells = [2764];%'all'; % 'all' or give a vector of vision cell ids


%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datarun{1}.names.rrs_params_path=['/Volumes/Lab/Projects/spikesorting/mvision/outputsSpectral/', run_opt.data_set, '/', run_opt.data_run, '/', 'data036', '.params'];

% datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data008-from-data009concate', '.params'];
% datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data022', '.sta'];
% datarun.names.rrs_neurons_path=['/Volumes/Lab/Projects/spikesorting/mvision/outputsSpectral/', run_opt.data_set, '/', run_opt.data_run, '/', 'data036', '.neurons'];

datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.neurons'];
% datarun.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', '08'];%run_opt.data_run(end-1:end)];
datarun.names.stimulus_path=['~/Desktop/2016-01-05-5/','s03','.txt'];


% datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.params'];
% datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.sta'];
%
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.neurons'];
% datarun.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', run_opt.data_run(end-1:end)];
opt=struct('verbose',1,'load_sta', 1,'load_params',0,'load_neurons',1,'load_obvius_sta_fits',true);

% datarun{1}=load_data(datarun{1},opt);
datarun=load_data(datarun,opt);

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
end


% get spikes per bin

% These triggers depend on the Labview trigger. 7.5 seconds prevents
% triggers from being missed on the large spatial scales
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

for x = 1:length(cells)
    fprintf('.');
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
            hh=find(h>=start & h<=stop); % Find spikes for each grating that passed over the spot
            psth_r=[psth_r; (h(hh))', repmat(length(sub_spacing)-i,[length(hh),1]);];
            
            
        end
%         rasterplot(psth_r(:,1)+ (psth_r(:,2)-1)*stop_time,interval/stop_time ,stop_time)
        psth_mat{z} = psth_r; % Put all the psths from all the temporal, spatial and direction combinations into a cell array
        %                if rgb(1) == 0.48
        %                    psth_mat{z} = [];
        %                end
        
        psth_r = [];
        
    end
    
    
    %     for i = 25:size(tr,2)
    %         figure
    %         plot(psth_mat{i}(:,1), psth_mat{i}(:,2), 'o')
    %     end
    
    
    % Organize the data into [direction, spatial, temporal, number of spikes in
    % 5 second trial]
    by_trial = zeros(size(tr,2),4);
    for i = 1:size(tr,2)-1
        by_trial(i,1) = datarun.stimulus.trials(i).DIRECTION;
        by_trial(i,2) = datarun.stimulus.trials(i).SPATIAL_PERIOD;
        by_trial(i,3) = datarun.stimulus.trials(i).TEMPORAL_PERIOD;
        by_trial(i,4) = size(psth_mat{i},1); % number of spikes in that 5 seconds regardless of timing
        
    end
    
    % Plot the data separated out by either temporal or spatial period
    fig = figure;
    %     set(fig,'PaperPositionMode','auto')
    set(fig, 'PaperSize', [8 12], 'PaperPosition', [0 0 8, 12]);
    set(fig, 'Position', [0 0 1200 800]);
    %     set(fig, 'PaperUnits', 'inches')
    %     set(fig, 'PaperPosition', [0 0 12 8]);
    
    D = unique(by_trial,'rows');
    placement = 1;
    for t = 1:length(temporal_tested)
        spat64 = find(by_trial(:,3) == temporal_tested(t));
        
        for s = 1:length(spatial_tested)
            spac = find(by_trial(:,2) == spatial_tested(s));
            C = intersect(spat64, spac);
            
            subplot(length(temporal_tested),length(spatial_tested),placement)
            angles = unique(by_trial(C,1));
            to_plot = [];
            error_u = [];
            error_l = [];
            
            for ang = 1:length(angles)
                indicies = find(by_trial(:,1) == angles(ang));
                just_one_angle = intersect(C, indicies);
                to_plot(ang,:) = [angles(ang) mean(by_trial(just_one_angle,4))];
                error_u(ang, :) = [angles(ang) mean(by_trial(just_one_angle,4))+std(by_trial(just_one_angle,4))];
                error_l(ang, :) = [angles(ang) mean(by_trial(just_one_angle,4)) - std(by_trial(just_one_angle,4))];
                
            end
            to_plot = [to_plot; to_plot(1,:)];
            error_u = [error_u; error_u(1,:)];
            error_l = [error_l; error_l(1,:)];
            error_l(error_l < 0) = 0;
            %                 h = polar(deg2rad(error_u(:,1)), error_u(:,2), 'r.-');
            hold on
            %                 g = polar(deg2rad(error_l(:,1)), error_l(:,2), 'r.-');
            
            j = polar(deg2rad(to_plot(:,1)), to_plot(:,2), 'k-');
            
            
%             polarwitherrorbar(deg2rad(to_plot(1,:)),to_plot(2,:),error_b)
            %                 h = polar(deg2rad(by_trial(C,1)), by_trial(C,4), 'b.');
            %                 set( findobj(h, 'Type', 'line'),'MarkerSize',15);
            %                 set( findobj(g, 'Type', 'line'),'MarkerSize',15);
            %                 set( findobj(j, 'Type', 'line'),'MarkerSize',15);
            
            hold on
                        k = polar(deg2rad(error_u(:,1)), error_u(:,2), 'r-');
                        l = polar(deg2rad(error_l(:,1)), error_l(:,2), 'r-');

%                             hHiddenText = findall(gca,'type','text');
%                             Angles = 0 : 45 : 315;
%                             hObjToDelete = zeros( length(Angles)-4, 1 );
%                             k = 0;
%                             for ang = Angles
%                                 k  =k+1;
%                                 hObj = findall(hHiddenText,'string',num2str(ang));
%                                 hObjToDelete(k) = hObj;
%             
%                             end
%                             delete( hObjToDelete(hObjToDelete~=0) );
%             
            
            % say your labels have the following strings..
            hHiddenText = findall(gca,'type','text');
            rho_labels = {'  10' '  20' '  30', '  40', '  50'};
            for r=1:length(rho_labels)
                delete(findall(hHiddenText, 'string', rho_labels{r}))
            end
             axis equal
%              axis square
            
            
            title({['Temporal = ', num2str(temporal_tested(t))]; ['Spatial = ', num2str(spatial_tested(s))]})
            placement = placement + 1;
        end
        clear spat64
        
    end
    
    %
    %    for t = 1:length(spatial_tested)
    %
    %     spat64 = find(by_trial(:,2) == spatial_tested(t));
    %     for i = spat64
    %         subplot(2,max(length(temporal_tested), length(spatial_tested)),t + max(length(temporal_tested), length(spatial_tested)))
    %
    %         polar(deg2rad(by_trial(i,1)), by_trial(i,4), 'k.')
    %         hold on
    %     end
    %     title(['Spatial = ', num2str(spatial_tested(t))])
    %     xlabel('Direction')
    %     ylabel('Spikes count during stimulus')
    %
    %     clear spat64
    %    end
    
    
    
    suptitle({run_opt.data_set;['Cell ', num2str(cell_to_run)]});
    
    
    if ~exist(filepath)
        mkdir(filepath)
    end
    
    % set(gcf, 'renderer', 'opengl')
    
    %     export_fig([filepath, 'Cell_',num2str(cell_to_run)], '-pdf')
    print(fig,'-dpdf',[filepath, 'Cell_',num2str(cell_to_run)]);
%     close(gcf);
end
