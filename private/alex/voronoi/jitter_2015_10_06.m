
datarun = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data008/data008');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-8-0.48-11111-160x160.xml');

disp('done inputs loading')

total_stix = 160*160;
a = permute(inputs, [2 1 3]);
tmp = a(:);
adj = 0;
for i=1:duration
    tmp([i*total_stix+1+adj,i*total_stix+2+adj]) = -1;
    adj=adj+2;
end
tmp(tmp==-1) = [];
new_inputs = reshape(tmp,total_stix,duration);

inputs = new_inputs;
clear tmp a new_inputs


sta_params.length = 6;
sta_params.offset = 0;
all_sta = cell(1,length(datarun.cell_ids));


spike_array = uint8(zeros(length(datarun.cell_ids), duration));
for cellID = 1:length(datarun.cell_ids)    
    spikes = datarun.spikes{cellID};
    spikes=ceil((spikes-datarun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes>duration) = [];
    spikes(spikes<sta_params.length-sta_params.offset) = [];    
    spike_array(cellID, spikes) = 1;
end

load('/Volumes/Analysis/2015-10-06-2/jitter/shifts')

% parts of STAs
disp('starting sta calculation')
parts = 50;
for kkk = 1:parts:length(datarun.cell_ids)
    
    spike_array_tmp = spike_array(kkk:kkk+parts-1,:);

    full_inputs = zeros(321,321,sta_params.length);
    for i=1:sta_params.length-1
        tmp = reshape(inputs(:,i), 160, 160);
        tmp = imresize(tmp,2, 'method', 'nearest');
        jitterX = shifts(1,i);
        jitterY = shifts(2,i);
        full_inputs(2+jitterX:321+jitterX,2+jitterY:321+jitterY,1+i) = tmp;
    end
    
    sta = zeros(321,321, sta_params.length,parts);
    tic
    for i=sta_params.length:duration
        tmp = reshape(inputs(:,i),160,160);
        tmp = imresize(tmp,2, 'method', 'nearest');
        jitterX = shifts(1,i);
        jitterY = shifts(2,i);
        
        full_inputs = circshift(full_inputs,-1,3);
        full_inputs(2+jitterX:321+jitterX,2+jitterY:321+jitterY,sta_params.length) = tmp;
        
        a = find(spike_array_tmp(:,i));
        sta(:,:,:,a) =  sta(:,:,:,a) + repmat(full_inputs, 1, 1, 1, length(a));
        if mod(i,10000)==0
            a = toc;
            disp(['processed ', int2str(i), ' frames in ', num2str(a), ' s'])
            tic;
            i
        end
    end

    for i = 1:parts
        sta(:,:,:,i) =  sta(:,:,:,i) / nnz(spike_array_tmp(i,1:duration));
    end
    save(['/Volumes/Analysis/2015-10-06-2/jitter/correct_jitter_sta_cells_',int2str(kkk),'.mat'], 'sta', '-v7.3');
    
end
% 
% 
% load('/Volumes/Analysis/2015-10-06-2/jitter/correct_jitter_sta_cells_1.mat')
% 
% for j=1:3
%     figure
%     for i=1:6
%         tmp = sta(:,:,i,j);
%         tmp = tmp/max(abs(tmp(:)))/2+0.5;
%         subplot(2,3,i)
%         imagesc(tmp);
%     end
% end
% 
% 
% for j=21:20
%     figure
%     for i=27
%         tmp = sta(:,:,:,i,j);
%         tmp = tmp/max(abs(tmp(:)))/2+0.5;
%         imagesc(tmp);
%     end
% end
% 
% 
