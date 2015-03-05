clear

chirpLength=22500;

%% get firing rates from 2 experiments
cd('S:\user\alexandra\scripts')
date='20120928';
path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    cnt=1;
    for i=10:12:82
        ampl(unit,cnt)=mean(common_res_fit(1,i:i+2));
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
    for i=1:5:35
        responses(cnt,1:chirpLength,unit)=mean(chirp(i:i+4,1:chirpLength));
        cnt=cnt+1;
    end

end

date='20121002';
path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    cnt=1;
    for i=22:12:94
        ampl(unit+24,cnt)=mean(common_res_fit(1,i:i+2));
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
    for i=6:5:40
        responses(cnt,1:chirpLength,unit+24)=mean(chirp(i:i+4,1:chirpLength));
        cnt=cnt+1;
    end
end

clear name common_res_fit chirp chirpPath date flips fits path2fit trial trialsChirp unit spike_info conv spikes

%% get parameters

noStim=200:2500;
%param 2-17 for fc
FC_beg=3000;
FC_step=500;
FC_end=11000;
% param 18:33 for ac
AC_beg=11670;
AC_step=500;
AC_end=19670;
% structure: rows - parameters,columns - nds, slices - units
for i=1:53
    for j=1:7
        baseLine=mean(responses(j,noStim,i));
        based=abs(responses(j,:,i)-baseLine);
        param(1,j,i)=sum(based(noStim))/length(noStim);
        cnt=2;
        for fc=FC_beg:FC_step:FC_end-500
            param(cnt,j,i)=sum(based(fc:fc+FC_step))/FC_step;
            cnt=cnt+1;
        end
        for ac=AC_beg:AC_step:AC_end-500
            param(cnt,j,i)=sum(based(ac:ac+AC_step))/AC_step;
            cnt=cnt+1;
        end
    end
end

clear based baseLine noStim lowFreqChirp midFreqChirp highFreqChirp lowAmplChirp midAmplChirp highAmplChirp

%% normalize amplitude
ampl=abs(ampl);
for i=1:53
    tmp=max(ampl(i,:));
    ampl(i,:)=ampl(i,1:7)/tmp;    
end

%% normalize to spontaneous modulation
spont=repmat(param(1,:,:),33,1);
% if spontaneous activity was 0, leave the response as it was
spont(spont==0)=1;
%get ratio to spontaneous modulation
param=param./spont;
clear spont

%% correlation for each unit between amplitude and parameter changes
k=zeros(53,33);
for j=1:53
    for i=1:33        
        k(j,i)=corr(ampl(j,:)',param(i,:,j)');
    end
end
nanmean(k)


figure
plot(1, (nanmean(k(:,1))),'.','markersize',20)
hold on
plot(2:17, (nanmean(k(:,2:17))),'linewidth',3)
plot(2:17, (nanmean(k(:,2:17))),'.','markersize',20)
plot(18:33, (nanmean(k(:,18:33))),'linewidth',3)
plot(18:33, (nanmean(k(:,18:33))),'.','markersize',20)
title('correlation between response modulation at different chirp intervals changes and filter amplitude changes across ND')
line([0 34],[0,0],'color','k')













%% FULL

noStim=200:2500;
window=500;
% structure: rows - parameters,columns - nds, slices - units
for i=1:53
    for j=1:7
        baseLine=mean(responses(j,noStim,i));
        based=abs(responses(j,:,i)-baseLine);
        cnt=1;
        for fc=1:window:19500
            param(cnt,j,i)=sum(based(fc:fc+window))/window;
            cnt=cnt+1;
        end
    end
end

clear based baseLine noStim lowFreqChirp midFreqChirp highFreqChirp lowAmplChirp midAmplChirp highAmplChirp

%% normalize amplitude
ampl=abs(ampl);
for i=1:53
    tmp=max(ampl(i,:));
    ampl(i,:)=ampl(i,1:7)/tmp;    
end

%% normalize to spontaneous modulation
spont=repmat(param(1,:,:),39,1);
% if spontaneous activity was 0, leave the response as it was
spont(spont==0)=1;
%get ratio to spontaneous modulation
param=param./spont;
clear spont

%% correlation for each unit between amplitude and parameter changes
k=zeros(53,39);
for j=1:53
    for i=1:39        
        k(j,i)=corr(ampl(j,:)',param(i,:,j)');
    end
end
nanmean(k)

figure
plot(1:39, (nanmean(k(:,1:39))),'linewidth',3)
hold on
plot(1:39, (nanmean(k(:,1:39))),'.','markersize',20)
line([0 39],[0,0],'color','k')
title('correlation between response modulation at different chirp intervals changes and filter amplitude changes across ND, FULL chirp response')






%% SLIDING WINDOW
clear param spont
noStim=200:2500;
window=500;
% structure: rows - parameters,columns - nds, slices - units
for i=1:53
    for j=1:7
        baseLine=mean(responses(j,noStim,i));        
        based=abs(responses(j,:,i)-baseLine);
% based=responses(j,:,i);
        spont(j,i)=sum(based(noStim));
        cnt=1;
        for fc=1:50:chirpLength-window
            param(cnt,j,i)=sum(based(fc:fc+window))/window;
            cnt=cnt+1;
        end        
        spont(1:cnt-1,j,i)=repmat(sum(based(noStim)),cnt-1,1);
    end
end

clear based baseLine noStim lowFreqChirp midFreqChirp highFreqChirp lowAmplChirp midAmplChirp highAmplChirp

%% normalize amplitude
ampl=abs(ampl);
for i=1:53
    tmp=max(ampl(i,:));
    ampl(i,:)=ampl(i,1:7)/tmp;    
end

%% normalize to spontaneous modulation
% if spontaneous activity was 0, leave the response as it was
spont(spont==0)=1;
%get ratio to spontaneous modulation
param=param./spont;
clear spont

%% correlation for each unit between amplitude and parameter changes
k=zeros(53,size(param,1));
for j=1:53
    for i=1:size(param,1)        
        k(j,i)=corr(ampl(j,:)',param(i,:,j)');
    end
end
nanmean(k)
h=k;
h(isnan(h))=1;

figure
plot(nanmean(k),'linewidth',3)
hold on
line([0 size(param,1)],[0,0],'color','k')
title('correlation between response modulation at different chirp intervals changes and filter amplitude changes across ND, FULL chirp response, sliding window 1ms')

figure
imagesc(cdata);
set(gca,'YDir','normal')
colormap(gray);
hold on
plot(std(h)+mean,'b','LineWidth',2);
plot(mean(h)+2,'r','LineWidth',2);
plot(min(h)+2,'g','LineWidth',2);
axis([0 chirpLength  1 5]);
line([0 chirpLength],[2 2],'color','k')
set(gca,'ytick',[1.8 2 2.3 2.7 3],'yticklabel',{'-0.2','0','0.3','0.7','1'})


figure
hold on
plot(std(h)./sqrt(53)+mean(h),'b','LineWidth',2);
plot(mean(h),'r','LineWidth',2);
plot(-std(h)./sqrt(53)+mean(h),'b','LineWidth',2);
axis([0 400  0 1]);

figure
plot(h');
figure
for i=1:7
    subplot(7,1,8-i)
    plot(mean(param(:,i,:),3))
    title(['ND',int2str(8-i),'   rod-only, n=53 (2 retinas)'])
end

param_rod=param;
save('C:\Documents and Settings\atikidzhi\My Documents\Dropbox\paper\2_lightAdaptation\chirp\rod_modul_param','param_rod')



rods=h;
save('C:\Documents and Settings\atikidzhi\My Documents\Dropbox\paper\2_lightAdaptation\chirp\rods_ampl_modul_corr','rods')