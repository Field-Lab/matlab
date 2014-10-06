% use EI to refine spike sorting

% PARAMETERS

% raw data path
raw_data_path = '/Volumes/War/Data/Gauthier/2007-09-18-2/data002/';

% spike window
spike_window = [0 30];
% sampling rate
sampling_rate = 20000;


% cell id to get template
%cell_id = 3080;
cell_id = 4066;

% electrodes
clear electrodes
switch 1
    case 1
        %electrodes{1} = [206 199 198]; % cell 3080
        electrodes{1} = [344 352 360 368 371 376 384 425 433 434 435 427 428 436 437 438 446 447 448 455 456]; % cell 4066
        electrodes{2} = [29 21 13 5 457 458 459 10 2 465 466 467 468 469 478 479 487]; % cell 4066
    case 2
        % choose electrodes within specified polygon
        xv = [98 98 968 968];yv = [0 177 777 0]; % cell 4066
        electrodes{1} = find(inpolygon(ep(:,1),ep(:,2),xv,yv));
    case 3
        electrodes{1} = 1:512;
        electrodes{2} = 1:512;
end



% LOAD TEMPLATE

%datarun = load_ei(datarun,cell_id);
ei = datarun.ei.eis{get_cell_indices(datarun,cell_id)};
clear template
for ee = 1:length(electrodes)
    %template{ee} = new_ei{ee}(electrodes{ee},:);
    template{ee} = ei(electrodes{ee},(spike_window(1):spike_window(2)) + datarun.ei.nlPoints + 1);
end



% PROJECT SPIKE WAVEFORMS ONTO TEMPLATE

% load data object
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);

% get spike times
spike_times = datarun.spikes{get_cell_indices(datarun,cell_id)};

% initialize projections
proj = zeros(length(spike_times),length(electrodes));

clear my_ei
for ee = 1:length(electrodes)
    my_ei{ee} = zeros(512,size(template{ee},2));
end
ei_count = zeros(ee,1);

% loop through spikes
for ss = 1:length(spike_times)
    % load voltage waveform on all electrodes
    all_waveforms = rawFile.getData(spike_times(ss)*sampling_rate + spike_window(1) + 1,...
        spike_window(2)-spike_window(1)+1);
    
    for ee = 1:length(electrodes)
        % pare to electrodes of interest
        waveforms = double(all_waveforms(:,electrodes{ee}+1)');

        % subtract mean
        waveforms = waveforms - repmat(mean(waveforms,2),1,size(waveforms,2));

        % project onto template
        proj(ss,ee) = sum(sum(waveforms .* template{ee}))/norm(reshape(waveforms,[],1));
        %proj(ss,ee) = norm(reshape(waveforms-template{ee},[],1));
            
        if mod(ss,100)==0 && ee == 1 && 0
            figure(4);clf;plot(template{ee}(:,:)');title('ei')
            figure(5);clf;plot(waveforms(:,:)');title('current waveform')
            figure(6);clf;plot(my_ei{ee}(electrodes{ee},:)'/length(spike_times));title('average spike')
            figure(7);clf;plot(reshape(template{ee},[],1),reshape(waveforms,[],1),'.');title('obs vs expected')
            hold on; plot([-50 50],[-50 50],'r');set(gca,'xlim',[-12 10],'ylim',[-150 100]);
            drawnow
        end
    end

    % add to ei 1
    if proj(ss,1)>20 && proj(ss,2)<0
        my_ei{1} = my_ei{1} + double(all_waveforms(:,2:513)') -...
            repmat(mean(double(all_waveforms(:,2:513)),1),size(all_waveforms,1),1)';
        ei_count(1) = ei_count(1) + 1;end
    % add to ei 2
    if proj(ss,1)<5 && proj(ss,2)>30
        my_ei{2} = my_ei{2} + double(all_waveforms(:,2:513)') - ...
            repmat(mean(double(all_waveforms(:,2:513)),1),size(all_waveforms,1),1)';
        ei_count(2) = ei_count(2) + 1;end

    % plot projections
    if mod(ss,100)==0
        figure(20);clf;
        plot(proj(:,1),proj(:,end),'.');
        %hist(proj(proj(:,1)~=0,2),100);drawnow
    end
end

% normalize EIs
my_ei{1} = my_ei{1}/ei_count(1);
my_ei{2} = my_ei{2}/ei_count(2);


