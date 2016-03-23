
datarun = load_data('/Volumes/Analysis/2015-10-06-2/d00-13-norefit/data008/data008');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-8-0.48-11111-160x160.xml');

disp('done inputs loading')

total_stix = 160*160;
new_inputs = zeros(total_stix,3,duration);
for j=1:3
    a = squeeze(inputs(:,j,:));
    a = reshape(a,160,160, duration);
    a = permute(a, [2 1 3]);
    tmp = a(:);
    adj = 0;
    for i=1:100000
        tmp([i*total_stix+1+adj,i*total_stix+2+adj]) = -1;
        adj=adj+2;
    end
    tmp(tmp==-1) = [];
    new_inputs(:,j,:) = reshape(tmp,total_stix,duration);
end
inputs = new_inputs;
clear tmp a new_inputs

k = unique(inputs(1:100));
inputs(inputs==k(2)) = 0.48;
inputs(inputs==k(1)) = -0.48;


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

    full_inputs = zeros(320+15,320+15,3,sta_params.length);
    for i=1:sta_params.length-1
        tmp = reshape(inputs(:,:,i), 160, 160, 3);
        tmp = imresize(tmp,2, 'method', 'nearest');
        jitterX = shifts(1,i);
        jitterY = shifts(2,i);
        full_inputs(9+jitterX:320+8+jitterX,9+jitterY:320+8+jitterY,:,1+i) = tmp;
    end
    
    sta = zeros(320+15,320+15, sta_params.length,parts);
    tic
    for i=sta_params.length:duration
        tmp = reshape(inputs(:,:,i),160,160, 3);
        tmp = imresize(tmp,2, 'method', 'nearest');
        jitterX = shifts(1,i);
        jitterY = shifts(2,i);
        
        full_inputs = circshift(full_inputs,-1,4);
        full_inputs(9+jitterX:320+8+jitterX,9+jitterY:320+8+jitterY,:,sta_params.length) = tmp;
        
        a = find(spike_array_tmp(:,i));
        sta(:,:,:,:,a) =  sta(:,:,:,:,a) + repmat(full_inputs, 1, 1, 1, 1, length(a));
        if mod(i,5000)==0
            a = toc;
            disp(['processed ', int2str(i), ' frames in ', num2str(a), ' s'])
            tic;
            i
        end
        if mod(i,40000)==0
            save(['/Volumes/Analysis/2015-10-06-2/jitter/correct_jitter_sta_',int2str(i),'_cells_',int2str(kkk),'.mat'], 'sta', '-v7.3');
            disp(['saved STA from ', int2str(i), ' frames for cells ', int2str(kkk), ' to ', int2str(kkk+parts-1)]);
        end
    end

    for i = 1:parts
        sta(:,:,:,:,i) =  sta(:,:,:,:,i) / nnz(spike_array_tmp(i,1:90000));
    end
    save(['/Volumes/Analysis/2015-10-06-2/jitter/correct_jitter_sta_complete_cells_',int2str(kkk),'.mat'], 'sta', '-v7.3');
    
end
% 
% 
% load('/Volumes/Analysis/2015-10-06-2/jitter/correct_jitter_sta_40000_cells_1.mat')
% 
% for j=1:3
%     figure
%     for i=1:30
%         tmp = sta(:,:,:,i,j);
%         tmp = tmp/max(abs(tmp(:)))/2+0.5;
%         subplot(5,6,i)
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
