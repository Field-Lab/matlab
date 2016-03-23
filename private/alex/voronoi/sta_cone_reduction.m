path2load = ['/Volumes/Analysis/2016-02-17-7/d01-13-norefit/data002/data002'];

datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);


path2load = ['/Volumes/Analysis/2016-02-17-7/d01-13-norefit/data003/data003'];

datarun1 = load_data(path2load);
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1,'load_sta','all','keep_java_sta',true);
datarun1 = load_neurons(datarun1);


for i=get_cell_indices(datarun, {4})
    sta = datarun.stas.stas{i}(:,:,:,4);
    sta = sta/max(abs(sta(:)));
    if -min(sta(:))>max(sta(:))
        sta = -sta;
    end
    sta = imresize(sta,2,'method', 'nearest');
    
    sta1 = datarun1.stas.stas{i}(:,:,:,4);
    sta1 = sta1/max(abs(sta1(:)));
    if -min(sta1(:))>max(sta1(:))
        sta1 = -sta1;
    end
    sta1 = imresize(sta1,3,'method', 'nearest');
    
    comb = zeros(600,800,3);
    comb(:,:,1) = sta;
    comb(:,1:795,2) = sta1;
    figure
    imagesc(comb)    
    
end


[inputs, refresh, duration] = get_wn_movie_ath(datarun1, 'BW-3-6-0.48-11111-265x200-60.35.xml');

inputs1 = inputs;
inputs = reshape(inputs, 200*265,[]);

cell_id = 75;%5854

spikes=ceil((datarun1.spikes{cell_id}-datarun1.triggers(1))*1000/(refresh)); % spikes in frames
    
sta_params.length = 6;
bin_nonl = 10; % number of bins for nonlinearity

spikes(spikes>size(inputs,2)) = []; 
spikes(spikes<sta_params.length) = [];
% get spike rate
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp
spike_rate(end+1:end+2) = 0;

% common STA to estimate bins for nonlinearity
spikes_tmp = spikes;
spikes_tmp(spikes_tmp<sta_params.length) = []; 
sta=zeros(200*265,sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    numel(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(inputs(:,spikes_tmp(ia) - sta_params.length + j),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;


figure
for i=1:6
    subplot(2,3,i)
    tmp = sta(:,i);
    imagesc(reshape(tmp,200,265));
end

frame = 2;
figure
tmp = reshape(sta(:,2), 200, 265);
imagesc(tmp);

tmp(37,88)
ind = find(sta(:,2)==tmp(37,88));

ind = find(sta(:,2)==tmp(32,85));

ind = find(sta(:,2)==tmp(35,88),1);

% ind = 1
figure
hold on
for j=1:length(ind)
    inp_cone = inputs(ind(j),:);
    inp_cone = [inp_cone(1:end-1); inp_cone(2:end)]';
    td = repmat(sta(ind(j),3:-1:2), 20733,1);
    inp_cone = inp_cone.*td;
    gen_signals = sum(inp_cone,2);
    
    a = unique(gen_signals);
    clear c
    for i=1:length(a)
        b = find(gen_signals==a(i));
        c(i)=mean(spike_rate(b+2));
    end    
    plot(a, c)
end
line([a(1), a(end)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
tt= std(spike_rate)/sqrt(length(spike_rate));
line([a(1), a(end)], [mean(spike_rate)+tt mean(spike_rate)+tt], 'color', 'k')




for j=1:length(ind)
    inp_cone = inputs(ind(j),:);
    td = repmat(sta(ind(j),2), 20734,1);
    inp_cone = inp_cone'.*td;
    gen_signals = inp_cone;
    
    a = unique(gen_signals);
    clear c
    for i=1:length(a)
        b = find(gen_signals==a(i));
        c(i)=mean(spike_rate(b+1));
    end    
    plot(a, c)    
end



for j=1:length(ind)
    inp_cone = inputs(ind(j),:);
    inp_cone = [inp_cone(1:end-2); inp_cone(2:end-1); inp_cone(3:end)]';
    td = repmat(sta(ind(j),3:-1:1), 20732,1);
    inp_cone = inp_cone.*td;
    gen_signals = sum(inp_cone,2);
    
    a = unique(gen_signals);
    clear c
    for i=1:length(a)
        b = find(gen_signals==a(i));
        c(i)=mean(spike_rate(b+2));
    end    
    plot(a, c)
end



b = [find(gen_signals==a(end)); find(gen_signals==a(1))];

spike_rate(b+2) = 0;
% get spike times from spike rate

spikes_tmp = find(spike_rate);

spikes_tmp(spikes_tmp<sta_params.length) = []; 
sta=zeros(200*265,sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    numel(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(inputs(:,spikes_tmp(ia) - sta_params.length + j),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;


figure
for i=1:6
    subplot(2,3,i)
    tmp = sta(:,i);
    imagesc(reshape(tmp,200,265));
end


figure
tmp = reshape(sta(:,2), 200, 265);
imagesc(tmp);



