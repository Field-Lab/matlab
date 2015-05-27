%%
datarun = load_data('/Volumes/Analysis/2015-04-14-2/data002_1_300_from_data000/data002_1_300_from_data000');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

wn_movie_name = 'RGB-20-1-0.48-11111-16x16.xml';

triggers=datarun.triggers;
mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', wn_movie_name];
[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,triggers, 1,2);
mvi=load_movie(mdf_file, triggers);


inputs=zeros(height,width,3,duration-1);

for j=1:duration
    F = round(mvi.getFrame(j-1).getBuffer);
    
    F = reshape(F,3,width,height);
    inputs(:,:,:,j) = permute(F,[3 2 1]);
end

inputs=(inputs*0.96)-0.48;


ncells = length(datarun.cell_ids);

offset = 0;
sta_length=20;

stas = cell(ncells,1);
nsp = zeros(ncells,1);
tic
parfor datarunID = 1:ncells
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/refresh); % spikes in frames

    spike_rate=zeros(duration,1);
    my_sta=zeros(16,16,3,sta_length);
    spikes(spikes<sta_length-offset)=[];
    nsp(datarunID) = length(spikes);
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        for j=1:sta_length
            my_sta(:,:,:,sta_length-j+1)=my_sta(:,:,:,sta_length-j+1)+sum(inputs(:,:,:,spikes(ia)-sta_length+j+offset),4);
        end
        spikes(ia)=[];
    end
%     my_sta=my_sta/nsp1(datarunID);
    stas{datarunID}=my_sta;
end
toc


%% 
datarun = load_data('/Volumes/Analysis/2015-04-14-2/data002_301_600_from_data000/data002_301_600_from_data000');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

wn_movie_name = 'RGB-20-1-0.48-11111-16x16.xml';


ncells = length(datarun.cell_ids);
t = [];
for datarunID = 1:ncells
    t = [t datarun.spikes{datarunID}(1)];
end
t = ceil((min(t)-datarun.triggers(1))*1000/refresh);

triggers=datarun.triggers;
mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', wn_movie_name];
[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,triggers, 1,2);
mvi=load_movie(mdf_file, triggers);

segment_movie_length = length((t-sta_length):duration);
inputs=zeros(height,width,3,segment_movie_length);

cnt = 1;
for j=(t-sta_length):duration
    F = round(mvi.getFrame(j-1).getBuffer);
    
    F = reshape(F,3,width,height);
    inputs(:,:,:,cnt) = permute(F,[3 2 1]);
    cnt = cnt+1;
end

inputs=(inputs*0.96)-0.48;


offset = 0;
sta_length=20;

stas1 = cell(ncells,1);
nsp1 = zeros(ncells,1);
tic
parfor datarunID = 1:ncells
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/refresh); % spikes in frames
    spikes = spikes-(t-sta_length); % get same starting point as movie inputs

    spike_rate=zeros(duration,1);
    my_sta=zeros(16,16,3,sta_length);
    spikes(spikes<sta_length-offset)=[];
    nsp1(datarunID) = length(spikes);
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        for j=1:sta_length
            my_sta(:,:,:,sta_length-j+1)=my_sta(:,:,:,sta_length-j+1)+sum(inputs(:,:,:,spikes(ia)-sta_length+j+offset),4);
        end
        spikes(ia)=[];
    end
%     my_sta=my_sta/nsp(datarunID);
    stas1{datarunID}=my_sta;
end
toc


%%
datarun = load_data('/Volumes/Analysis/2015-04-14-2/data002_601_900_from_data000/data002_601_900_from_data000');
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

wn_movie_name = 'RGB-20-1-0.48-11111-16x16.xml';

ncells = length(datarun.cell_ids);
t = [];
for datarunID = 1:ncells
    t = [t datarun.spikes{datarunID}(1)];
end
t = ceil((min(t)-datarun.triggers(1))*1000/refresh);

triggers=datarun.triggers;
mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', wn_movie_name];
[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,triggers, 1,2);
mvi=load_movie(mdf_file, triggers);

segment_movie_length = length((t-sta_length):duration);
inputs=zeros(height,width,3,segment_movie_length);

cnt = 1;
for j=(t-sta_length):duration
    F = round(mvi.getFrame(j-1).getBuffer);
    
    F = reshape(F,3,width,height);
    inputs(:,:,:,cnt) = permute(F,[3 2 1]);
    cnt = cnt+1;
end

inputs=(inputs*0.96)-0.48;


offset = 0;
sta_length=20;

stas2 = cell(ncells,1);
nsp2 = zeros(ncells,1);
tic
parfor datarunID = 1:ncells
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/refresh); % spikes in frames
    spikes = spikes-(t-sta_length); % get same starting point as movie inputs

    spike_rate=zeros(duration,1);
    my_sta=zeros(16,16,3,sta_length);
    spikes(spikes<sta_length-offset)=[];
    nsp2(datarunID) = length(spikes);
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        for j=1:sta_length
            my_sta(:,:,:,sta_length-j+1)=my_sta(:,:,:,sta_length-j+1)+sum(inputs(:,:,:,spikes(ia)-sta_length+j+offset),4);
        end
        spikes(ia)=[];
    end
%     my_sta=my_sta/nsp(datarunID);
    stas2{datarunID}=my_sta;
end
toc


%%

figure
for datarunID = 1:16%ncells
%     total_nsp = nsp(datarunID) + nsp1(datarunID) + nsp2(datarunID);
%     tt = (stas{datarunID}+stas1{datarunID}+stas2{datarunID})/total_nsp;
%     
    total_nsp = nsp2(datarunID);
    tt = stas2{datarunID}/total_nsp;

    sta = tt(:,:,:,6);
    tt = 0.5/max(abs(sta(:)));
    
    sta = sta*tt+0.5;
    subplot(4,4,datarunID)
    imagesc(sta)
    title(['Cell ', int2str(datarun.cell_ids(datarunID)), ', n spikes = ', int2str(total_nsp)])
    
end



figure
datarunID = 5;

%     total_nsp = nsp(datarunID) + nsp1(datarunID) + nsp2(datarunID);
%     tt = (stas{datarunID}+stas1{datarunID}+stas2{datarunID})/total_nsp;
%     
total_nsp = nsp1(datarunID);
sta = stas1{datarunID}/total_nsp;
tt = 0.5/max(abs(sta(:)));
sta = sta*tt+0.5;
for i=1:9
    tt = sta(:,:,:,i);
    subplot(3,3,i)
    imagesc(tt)
    title(['Cell ', int2str(datarun.cell_ids(datarunID)), ', n spikes = ', int2str(total_nsp)])
    
end



