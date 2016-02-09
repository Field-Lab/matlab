stimulus = zeros(50);

 % take STAs from a real cell (with good cones)? Gaussian around? 
 % For single cones. Adjust weights for Gauss or equal, or leave alone
 linear_filter = x;

% some function. Ignore subunits. Pretend Bipolar cell sums up all inputs 
% linearly, and nonlinearity between BC/RGC and own RGC is just one nonlinearity.
rgc_nonlin = y;

% add some poisson process to generate spikes.
% Then do reverse correlation to get STA. COmpare to initial linear
% filters. 
% add noise next step and repeat.



datarun = load_data('/Volumes/Analysis/2010-09-24-1/data002/data002');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun, 'load_sta', [1351]);
datarun = load_neurons(datarun);

datarunID = 35;
sta = squeeze(datarun.stas.stas{datarunID});

figure
colormap gray
imagesc(sta(:,:,5))

sta_tmp = sta(:,:,5);
flag = 1;
i = 1;
clear cones
while flag
    k = min(sta_tmp(:));
    [c,r] = find(sta_tmp==k);
    sta_tmp(c-1:c+1,r-1:r+1) = 0;
    cones(i,1) = r;
    cones(i,2) = c;
    i = i+1;
    if k>-0.020
        flag = 0;
    end
end

figure
colormap gray
imagesc(sta(:,:,5))
hold on
plot(cones(:,1),cones(:,2),'xg')

figure
hold on
clear mySTA
for i=1:length(cones)
    mySTA(i,:) = squeeze(sta(cones(i,2),cones(i,1),:));
    plot( mySTA(i,:))
end

myCones = cones;
myCones(:,1) = myCones(:,1) - min(cones(:,1))+5;
myCones(:,2) = myCones(:,2) - min(cones(:,2))+5;

t = sub2ind([40,40], myCones(:,1),myCones(:,2));
fullSTA = rand(1600,6)/200;
fullSTA(t,:) = mySTA;

t1 = sub2ind([320,320], cones(:,2),cones(:,1));
[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-1-6-0.48-11111-320x320.xml');

sta_params.length = 5;
sta_params.offset = 0;
spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];
fraction = 0.9;

[unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(t1,:), spikes, fraction, sta_params);

figure
plot(unbiased_sta')

figure
for i=1:25
    subplot(5,5,i)
    plot(gensig_bins(i,1:10), nonlinearity(i,:))
end





[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-1-6-0.48-11111-40x40.xml');

filtered_inputs=zeros(size(inputs,1),size(inputs,2)-6+1);  
for i=1:size(inputs,1)
    filtered_inputs(i,:)=conv(inputs(i,:),fullSTA(i,:),'valid');
end

gs = sum(filtered_inputs);
fr = zeros(size(gs));
for i=1:10
    a = find(gs>=gensig_bins(1,i) & gs<gensig_bins(1,i+1));
    fr(a) = nonlinearity(i);
end
dt = 1;
nBins = length(fr);
spikeMat = rand(1, nBins) < fr*dt;


spikes = find(spikeMat);
spikes(spikes<sta_params.length-sta_params.offset) = [];

spikes_tmp = spikes;
sta=zeros(size(inputs,1),sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    length(ia)
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(inputs(:,spikes_tmp(ia) - sta_params.length + j + sta_params.offset),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;

figure
for i=1:5
    subplot(2,3,i)
    colormap gray
    imagesc(reshape(sta(:,i), 40,40))
end


%% coarse
% data038 Binary BW 8-8-0.48-11111 ndf 3.9 40x40 1800 s
datarun = load_data('/Volumes/Analysis/2010-09-24-1/data038/data038');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun, 'load_sta', [2976]);
datarun = load_neurons(datarun);

datarunID = 79;
sta = squeeze(datarun.stas.stas{datarunID});

figure
colormap gray
imagesc(sta(:,:,26))

[c,r] = find(sta(:,:,26)<-0.03);
cones = [ r c];
t1 = sub2ind([40,40], cones(:,2),cones(:,1));

figure
colormap gray
imagesc(sta(:,:,26))
hold on
plot(cones(:,1),cones(:,2),'xg')

n_cones = length(cones);

[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-8-8-0.48-11111-40x40.xml');

sta_params.length = 20;
sta_params.offset = 0;
fraction = 0.9;

spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];

% spike rate
spikes(spikes>size(inputs,2)) = []; 
spikes(spikes<sta_params.length) = [];
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp

% [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(t1,:), spikes, fraction, sta_params);

spikes_tmp = spikes;
spikes_tmp(spikes_tmp<sta_params.length) = []; 
sta=zeros(n_cones,sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(inputs(t1,spikes_tmp(ia) - sta_params.length + j + sta_params.offset),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;

figure
plot(sta')

bin_nonl = 10;
gensig_bins = zeros(n_cones, bin_nonl+1);
tmp_spike_rate = spike_rate(sta_params.length:end);
for current_cone=1:n_cones
    filt_inputs=conv(inputs(t1(current_cone),:), sta(current_cone,:),'valid');
    
    bin_size = ceil(size(filt_inputs,2)/bin_nonl);
    tmp = sort(filt_inputs);    
    gensig_bins(current_cone,:) = [tmp(1:bin_size:end) tmp(end)*1.001];
    
    my_nl=zeros(bin_nonl,1);
    % mean firing rate in each bin of generator signal
    for k=1:bin_nonl
        my_nl(k)=mean(tmp_spike_rate(filt_inputs>=gensig_bins(current_cone,k) & filt_inputs<gensig_bins(current_cone,k+1)));
    end
    subplot(5,5,current_cone)
    hold on
    plot(gensig_bins(current_cone,1:end-1),my_nl);
    
end


bin_nonl = 10;
tmp_spike_rate = spike_rate(sta_params.length:end);
suminputs = 0;
for current_cone=1:n_cones
    filt_inputs=conv(inputs(t1(current_cone),:), sta(current_cone,:),'valid');
    suminputs = suminputs+filt_inputs;
end
bin_size = ceil(size(filt_inputs,2)/bin_nonl);
tmp = sort(suminputs);
gensig_bins= [tmp(1:bin_size:end) tmp(end)*1.001];
    
my_nl=zeros(bin_nonl,1);
for k=1:bin_nonl
    my_nl(k)=mean(tmp_spike_rate(suminputs>=gensig_bins(k) & suminputs<gensig_bins(k+1)));
end

plot(gensig_bins(1:end-1),my_nl);
    
% now use nonlinearity and STAs to predict the response


filtered_inputs=zeros(size(inputs,1),size(inputs,2)-sta_params.length+1);  
for i=1:size(t1,1)
    filtered_inputs(t1(i),:)=conv(inputs(t1(i),:),sta(i,:),'valid');
end

gs = sum(filtered_inputs);
fr = zeros(size(gs));
for i=1:10
    a = find(gs>=gensig_bins(i) & gs<gensig_bins(i+1));
    fr(a) = my_nl(i);
end
dt = 1;
spikeMat = rand(3, length(fr)) < repmat(fr,3,1)*dt;

% spike rate
spike_rate2 = 0;
for i=1:3
    spikes = find(spikeMat(i,:));
    spikes(spikes>size(inputs,2)) = [];
    spikes(spikes<sta_params.length) = [];
    spikes_tmp = spikes;
    spike_rate1=zeros(size(inputs,2),1);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate1(spikes_tmp(ia))=spike_rate1(spikes_tmp(ia))+1;
        spikes_tmp(ia)=[];
    end
    clear spikes_tmp
    spike_rate2 = spike_rate2+spike_rate1;
end

tmp = spike_rate2;
p1 = reshape(tmp(1:27000),10,2700);
p1 = mean(p1);
figure
plot(p1)

figure
plot(spike_rate)
hold on
plot(spike_rate1(15:end))

tmp = spike_rate;
p = reshape(tmp(1:27000),10,2700);
p = mean(p);

figure
plot(p)
hold on
plot(p1/2)

figure
plot(xcorr(p',p1', 50, 'coeff'))


spikes(spikes<sta_params.length-sta_params.offset) = [];

spikes_tmp = spikes;
sta=zeros(size(inputs,1),sta_params.length);
nspikes = numel(spikes_tmp);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    length(ia)
    for j=1:sta_params.length
        sta(:,sta_params.length-j+1) = sta(:,sta_params.length-j+1)...
            + sum(inputs(:,spikes_tmp(ia) - sta_params.length + j + sta_params.offset),2);
    end
    spikes_tmp(ia)=[];
end
sta = sta/nspikes;

figure
for i=1:20
    subplot(4,5,i)
    colormap gray
    imagesc(reshape(sta(:,i), 40,40))
end