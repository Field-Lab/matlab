clear

%% get firing rates from CONE experiment
date='20120920';
path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    cnt=1;
    for i=1:12:72
        ampl(unit,cnt)=mean(common_res_fit(1,i:i+2));
        lat(unit,cnt)=mean(common_res_fit(2,i:i+2));
        cnt=cnt+1;
    end
    
    name=fits(unit).name(1:end-15);
%     CHIRP
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',name,'___spike_info_chirp'])
    chirp=zeros(trialsChirp,25000);
    for trial=1:trialsChirp
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;
    end
    cnt=1;
    for i=1:5:30
        responses(cnt,1:22500,unit)=mean(chirp(i:i+4,1:22500));
        cnt=cnt+1;
    end
end


date='20120920_1';
path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    cnt=1;
    for i=1:12:60
        ampl(unit+38,cnt)=NaN;
        lat(unit+38,cnt)=NaN;
    end
    cnt=2;
    for i=1:12:60
        ampl(unit+38,cnt)=mean(common_res_fit(1,i:i+2));
        lat(unit+38,cnt)=mean(common_res_fit(2,i:i+2));
        cnt=cnt+1;
    end
    
    name=fits(unit).name(1:end-15);
    % CHIRP
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',name,'___spike_info_chirp'])
    chirp=zeros(trialsChirp,25000);
    for trial=1:trialsChirp
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;
    end
    cnt=1;
    for i=1:5:25
        responses(cnt,1:22500,unit+38)=NaN;
    end
    cnt=2;
    for i=1:5:25
        responses(cnt,1:22500,unit+38)=mean(chirp(i:i+4,1:22500));
        cnt=cnt+1;
    end
end



date='20120921';
path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    cnt=1;
    for i=1:12:60
        ampl(unit+63,cnt)=NaN;
        lat(unit+63,cnt)=NaN;
    end
    cnt=2;
    for i=1:12:60
        ampl(unit+63,cnt)=mean(common_res_fit(1,i:i+2));
        lat(unit+63,cnt)=mean(common_res_fit(2,i:i+2));
        cnt=cnt+1;
    end
    
    name=fits(unit).name(1:end-15);
    % CHIRP
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',name,'___spike_info_chirp'])
    chirp=zeros(trialsChirp,25000);
    for trial=1:trialsChirp
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;
    end
    cnt=1;
    for i=1:5:25
        responses(cnt,1:22500,unit+63)=NaN;
    end
    cnt=2;
    for i=1:5:25
        responses(cnt,1:22500,unit+63)=mean(chirp(i:i+4,1:22500));
        cnt=cnt+1;
    end
end


date='20120925';
path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    cnt=1;
    for i=1:12:60
        ampl(unit+72,cnt)=mean(common_res_fit(1,i:i+2));
        lat(unit+72,cnt)=mean(common_res_fit(2,i:i+2));
        cnt=cnt+1;
    end
    cnt=6;
    for i=1:12:60
        ampl(unit+72,cnt)=NaN;
        lat(unit+72,cnt)=NaN;
    end
    name=fits(unit).name(1:end-15);
    % CHIRP
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',name,'___spike_info_chirp'])
    chirp=zeros(trialsChirp,25000);
    for trial=1:trialsChirp
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;
    end
    cnt=1;
    for i=1:5:25
        responses(cnt,1:22500,unit+72)=mean(chirp(i:i+4,1:22500));
        cnt=cnt+1;
    end
    cnt=6;
    for i=1:5:25
        responses(cnt,1:22500,unit+72)=NaN;
    end
end



clear name common_res_fit chirp chirpPath date flips fits path2fit trial trialsChirp unit spike_info conv spikes

%% SLIDING WINDOW
clear param spont
noStim=200:2500;
window=500;
% structure: rows - parameters,columns - nds, slices - units
for i=1:91
    for j=1:6
        baseLine=mean(responses(j,noStim,i));        
        based=abs(responses(j,:,i)-baseLine);
% based=responses(j,:,i);
        spont(j,i)=sum(based(noStim));
        cnt=1;
        for fc=1:50:22500-window
            param(cnt,j,i)=sum(based(fc:fc+window))/window;
            cnt=cnt+1;
        end        
        spont(1:cnt-1,j,i)=repmat(sum(based(noStim)),cnt-1,1);
    end
end

clear based baseLine noStim lowFreqChirp midFreqChirp highFreqChirp lowAmplChirp midAmplChirp highAmplChirp

%% normalize amplitude
ampl=abs(ampl);
for i=1:91
    tmp=max(ampl(i,:));
    ampl(i,:)=ampl(i,1:6)/tmp;    
end


%% correlation for each unit between amplitude and parameter changes
k=zeros(91,size(param,1));
for j=1:91
    for i=1:size(param,1)        
        k(j,i)=corr(ampl(j,:)',param(i,:,j)');
%         k(j,i)=corr(lat(j,:)',param(i,:,j)');
    end
end
nanmean(k)
h=k;
h(isnan(h))=1;

figure
plot(nanmean(k),'linewidth',3)


figure
hold on
plot(std(h)./sqrt(57)+mean(h),'b','LineWidth',2);
plot(mean(h),'r','LineWidth',2);
plot(-std(h)./sqrt(57)+mean(h),'b','LineWidth',2);
% axis([0 400  0 1]);

figure
plot(h');

figure
for i=1:5
    subplot(5,1,6-i)
    plot(mean(param(:,i,:),3))
    title(['ND',int2str(6-i),'   cone-only, n=72 (3 retinas)'])
end
param_cone=param;
save('C:\Documents and Settings\atikidzhi\My Documents\Dropbox\paper\2_lightAdaptation\chirp\cone_modul_param','param_cone')



cone=h;
save('C:\Documents and Settings\atikidzhi\My Documents\Dropbox\paper\2_lightAdaptation\chirp\cone_ampl_modul_corr','cone')
