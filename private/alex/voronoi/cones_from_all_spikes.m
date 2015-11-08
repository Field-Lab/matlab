

local_path = '/Volumes/Analysis/';


datarun = load_data([local_path, '2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11']);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

tic
[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-6-0.48-11111-300x300-60.35.xml');
toc

% get spikes in frames
spikes = cell2mat(datarun.spikes);
spikes = round(sort(spikes)*1000/refresh);
spikes(spikes<2) = [];

spikes_tmp = spikes;
spike_rate=zeros(1,max(spikes));
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp
figure
plot(spike_rate)

inputs(inputs<0) = 0;
inputs(inputs>0) = 1;
% inputs transitions
tic
a = diff(inputs,1,2);
toc
a = a(:,1:18000);

meanpos = zeros(1,90000);
meanneg = meanpos;
for ind = 1:90000
    pos=a(ind,:)>0; % pos transitions
    neg = a(ind,:)<0; % neg transitions
    meanpos(ind)=mean(spike_rate(pos))/sum(pos);
    meanneg(ind)=mean(spike_rate(neg))/sum(neg);
end
[k, ic] = sort(meanneg);
figure
plot(meanneg(ic))
hold on
plot(meanpos(ic))


robust_mean(meanpos)-3*robust_std(meanpos)
robust_mean(meanneg)-3*robust_std(meanneg)

tmp=meanneg+meanpos;

figure
colormap gray
imagesc(reshape(tmp,300,300))



spikes = datarun.spikes{11};
spikes = round(sort(spikes)*1000/refresh);
spikes(spikes<2) = [];

spikes_tmp = spikes;
spike_rate_tmp=zeros(1,max(spikes));
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate_tmp(spikes_tmp(ia))=spike_rate_tmp(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp
figure
plot(spike_rate_tmp)


meanpos_tmp = zeros(1,90000);
meanneg_tmp = meanpos_tmp;
for ind = 1:90000
    pos=a(ind,:)>0; % pos transitions
    neg = a(ind,:)<0; % neg transitions
    lal = [0 spike_rate_tmp];
    meanpos_tmp(ind)=mean(lal(pos))/sum(pos);
    meanneg_tmp(ind)=mean(lal(neg))/sum(neg);
end


tmp=meanneg+meanpos;

figure
colormap gray
imagesc(reshape(tmp,300,300))
