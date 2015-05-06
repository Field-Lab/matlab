% ds_cells
clear
run_opt.data_set = '2015-04-09-1';
run_opt.data_run = 3; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
if strcmp(run_opt.data_set, '2015-04-09-1')
    datarun{1}.names.rrs_params_path='/Volumes/Analysis/2015-04-09-1/data015-mapped-data011/data015-mapped-data011.params';
    datarun{2}.names.rrs_neurons_path=('/Volumes/Analysis/2015-04-09-1/data015-mapped-data011/data015-mapped-data011.neurons');
    datarun{2}.names.stimulus_path=('/Volumes/Data/2015-04-09-1/Visual/s15');
end
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);
datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0, 'trigger_iti_thr', 0.0006); % manually set threshold until got right number of triggers
% Something went wrong with the first two stimuli
% datarun{2}.stimulus.triggers= [datarun{2}.stimulus.triggers(1:33) 187.0086, datarun{2}.stimulus.triggers(34:70) 396.0097 datarun{2}.stimulus.triggers(71:72) 412.5063 datarun{2}.stimulus.triggers(73:72); % Remove the first trigger
% datarun{2}.stimulus.trials = datarun{2}.stimulus.trials(3:end); % Remove the first two trials


cell_to_run = 889;
        cell_indices1=get_cell_indices(datarun{1}, [cell_to_run]);
        cell_indices2=get_cell_indices(datarun{2},[cell_to_run]);

% get spikes per bin

     tr=datarun{2}.stimulus.triggers; % all start triggers
        run_opt.config_num = 1;

get_stimulus_triggers(datarun)
       
n = datarun{2}.spikes{cell_indices2(1)}';
psth_r = [];
for z = 1:size(tr,2)
        ind = find(tr(z) == datarun{2}.triggers);
        spatial = datarun{2}.stimulus.trials(z).SPATIAL_PERIOD;
        temporal = datarun{2}.stimulus.trials(z).TEMPORAL_PERIOD;
        direction = datarun{2}.stimulus.trials(z).DIRECTION;
        interval = 5;
        num_triggers = interval/(temporal/120);
        sub_spacing= datarun{2}.triggers(ind:ind+num_triggers) - datarun{2}.triggers(ind);
        start = 0;
        stop = temporal/120;

        

        
for i=1:length(sub_spacing)-1,
        h=n-tr(z)-sub_spacing(i);  
        stop = sub_spacing(i+1) - sub_spacing(i);
        hh=find(h>=start & h<=stop);
        psth_r=[psth_r; (h(hh)*1000)', repmat(length(sub_spacing)-i,[length(hh),1]);];
        
        
end
psth_mat{z} = psth_r;
psth_r = [];

end


k=1; kmin=1; kmax=size(psth_mat,2); hk=loop_slider(k,kmin,kmax);

% while k
%     if ~ishandle(hk)
%         break % script breaks until figure is closed
%     end
%         
%         k=round(get(hk,'Value'));
%         
%         plot(psth_mat{k}(:,1), psth_mat{k}(:,2), '.k')
%         title({['Cell: ', num2str(cell_to_run)]; ['Spatial Period: ', num2str(datarun{2}.stimulus.trials(k).SPATIAL_PERIOD)]; ['Temporal Period: ', num2str(datarun{2}.stimulus.trials(k).TEMPORAL_PERIOD)]; ['Direction: ', num2str(datarun{2}.stimulus.trials(k).DIRECTION)]})
% drawnow
% hold off
% uiwait;
% end

by_trial = zeros(size(tr,2),4);
for i = 1:size(tr,2)
    by_trial(i,1) = datarun{2}.stimulus.trials(i).DIRECTION;
        by_trial(i,2) = datarun{2}.stimulus.trials(i).SPATIAL_PERIOD;
        by_trial(i,3) = datarun{2}.stimulus.trials(i).TEMPORAL_PERIOD;
        by_trial(i,4) = size(psth_mat{i},1);

end

spat64 = find(by_trial(:,3) == 60);
figure;
for i = spat64
    plot(by_trial(i,1), by_trial(i,4), 'k.')
    hold on
end

spat64 = find(by_trial(:,3) == 30);
figure;
for i = spat64
    plot(by_trial(i,1), by_trial(i,4), 'k.')
    hold on
end

spat64 = find(by_trial(:,3) == 15);
figure;
for i = spat64
    plot(by_trial(i,1), by_trial(i,4), 'k.')
    hold on
end

spat64 = find(by_trial(:,2) == 64);
figure;
for i = spat64
    plot(by_trial(i,1), by_trial(i,4), 'k.')
    hold on
end

spat64 = find(by_trial(:,2) == 32);
figure;
for i = spat64
    plot(by_trial(i,1), by_trial(i,4), 'k.')
    hold on
end

spat64 = find(by_trial(:,2) == 16);
figure;
for i = spat64
    plot(by_trial(i,1), by_trial(i,4), 'k.')
    hold on
end

    