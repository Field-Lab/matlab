
datarun = load_data('/Volumes/Analysis/2016-02-17-6/data026/data026');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

[inputs, refresh, duration] = get_wn_movie_ath_rgb(datarun, 'RGB-16-2-0.48-22222-119.5.xml');

disp('done inputs loading')

new_inputs = zeros(800,3,107301);
for j=1:3
    a = squeeze(inputs(:,j,:));
    tmp = a(:);
    adj = 0;
    for i=1:100000
        tmp(i*800+1+adj:i*800+2+adj) = -1;
        adj=adj+2;
    end
    tmp(tmp==-1) = [];
    new_inputs(:,j,:) = reshape(tmp,800,107301);
end
inputs = new_inputs;
clear tmp a new_inputs

k = unique(inputs(1:100));
inputs(inputs==k(2)) = 0.48;
inputs(inputs==k(1)) = -0.48;


sta_params.length = 30;
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

load('/Volumes/Analysis/2016-02-17-6/jitter/shifts')

% parts of 20 STAs
disp('starting sta calculation')

for kkk = 1:20:length(datarun.cell_ids)
    
    spike_array_tmp = spike_array(kkk:kkk+49,:);
    
    full_inputs = zeros(320+15,640+15,3,sta_params.length);
    for i=1:sta_params.length-1
        tmp = reshape(inputs(:,:,i), 20,40, 3);
        tmp = imresize(tmp,16, 'method', 'nearest');
        jitterX = shifts(1,i);
        jitterY = shifts(2,i);
        full_inputs(9+jitterX:320+8+jitterX,9+jitterY:640+8+jitterY,:,1+i) = tmp;
    end
    
    sta = zeros(320+15,640+15, 3, sta_params.length,20);
    for i=sta_params.length:90000
        tmp = reshape(inputs(:,:,i),20,40, 3);
        tmp = imresize(tmp,16, 'method', 'nearest');
        jitterX = shifts(1,i);
        jitterY = shifts(2,i);
        
        full_inputs = circshift(full_inputs,-1,29);
        full_inputs(9+jitterX:320+8+jitterX,9+jitterY:640+8+jitterY,:,sta_params.length) = tmp;
        
        a = find(spike_array_tmp(:,i));
        sta(:,:,:,:,a) =  sta(:,:,:,:,a) + repmat(full_inputs, 1, 1, 1, 1, length(a));
        
        if mod(i,500) == 0
            kkk
            i
        end
        if mod(i,500)==0
            save(['/Volumes/Analysis/2016-02-17-6/jitter/correct_jitter_sta_',int2str(i),'_cells_',int2str(kkk),'.mat'], 'sta', '-v7.3');
        end
    end
    
    for i = 1:20
        sta(:,:,:,:,i) =  sta(:,:,:,:,i) / nnz(spike_array_tmp(i,1:90000));
    end
    save(['/Volumes/Analysis/2016-02-17-6/jitter/correct_jitter_sta_complete_cells_',int2str(kkk),'.mat'], 'sta', '-v7.3');
    
end
