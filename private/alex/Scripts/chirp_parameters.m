clear

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
        responses(cnt,1:20000,unit)=mean(chirp(i:i+4,1:20000));
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
        responses(cnt,1:20000,unit+24)=mean(chirp(i:i+4,1:20000));
        cnt=cnt+1;
    end
end

clear name common_res_fit chirp chirpPath date flips fits path2fit trial trialsChirp unit spike_info conv spikes

%% get parameters

noStim=200:2500;
lowFreqChirp=3030:4935;
% midFreqChirp=4935:7953;
% midFreqChirp=4935:9720;
midFreqChirp=4935:7953;
highFreqChirp=7953:10970;
lowAmplChirp=11670:12600;
midAmplChirp=12600:16100;
highAmplChirp=16100:19600;

m=['spo';'LFC';'MFC';'HFC';'LAC';'MAC';'HAC';'tot']; % parameters order
% structure: rows - parameters,columns - nds, slices - units
for i=1:53
    for j=1:7
        baseLine=mean(responses(j,noStim,i));
        based=abs(responses(j,:,i)-baseLine);
        param(1,j,i)=sum(based(noStim))/length(noStim);
        param(2,j,i)=sum(based(lowFreqChirp))/length(lowFreqChirp);
        param(3,j,i)=sum(based(midFreqChirp))/length(midFreqChirp);
        param(4,j,i)=sum(based(highFreqChirp))/length(highFreqChirp);
        param(5,j,i)=sum(based(lowAmplChirp))/length(lowAmplChirp);
        param(6,j,i)=sum(based(midAmplChirp))/length(midAmplChirp);
        param(7,j,i)=sum(based(highAmplChirp))/length(highAmplChirp);
        param(8,j,i)=sum(responses(j,:))/length(responses);
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
spont=repmat(param(1,:,:),8,1);
% if spontaneous activity was 0, leave the response as it was
spont(spont==0)=1;
%get ratio to spontaneous modulation
param=param./spont;
clear spont

%% correlation for each unit between amplitude and parameter changes
k=zeros(53,8);
for j=1:53
    for i=1:8        
        k(j,i)=corr(ampl(j,:)',param(i,:,j)');
    end
end
nanmean(k)
% without normalizing to spontaneous
% 0.5011    0.7735    0.7883    0.8751    0.7667    0.7883    0.6997    0.8394
% with normalizing to spontaneous
% 0.3225    0.3188    0.3524    0.7182    0.6707    0.5118    0.2753    0.3774
for nd=1:7
    for i=1:8
        subplot(2,4,i)
        hold on
        tmp=param(i,nd,:)<50;
        tmp=tmp(:);
        plot(ampl(:,nd),reshape(param(i,nd,:),1,53),'.')
        plot(nanmean(ampl(:,nd)),nanmean(param(i,nd,tmp)),'k.','Markersize',20)
        
        tmp1=param(i,4,:)<50;
        tmp1=tmp1(:);
        plot(ampl(:,4),reshape(param(i,4,:),1,53),'r.')
        plot(nanmean(ampl(:,4)),nanmean(param(i,4,tmp1)),'k.','Markersize',20)
        axis([-0.5 1.5 0 50])
        title(m(i,:))
        
        tmp=tmp&tmp1;
        
        [h(nd,i),p(nd,i)]=ttest(reshape(param(i,2,tmp),1,sum(tmp))',reshape(param(i,4,tmp),1,sum(tmp))');
        %     [~,t(i)]=ttest(reshape(param(i,2,:),1,53)',reshape(param(i,4,:),1,53)');
        
    end
end
figure
plot(p','*')

l=p;
l=[(7:-1:1)',l]
l=[0:8; l]



figure
nds='7654321'
clear t
for nd=1:7
    subplot(2,4,nd)
    for i=1:8
        tmp=param(i,nd,:)<50;
        tmp=tmp(:);
        t(i,1)=nanmean(param(i,nd,tmp));
        t(i,2)=nanstd(param(i,nd,tmp))/sqrt(sum(~isnan(param(i,nd,tmp))));
    end
    errorbar(t(:,1),t(:,2))
    title(['ND',nds(nd)])    
    set(gca,'XTick',1:8,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    axis([0 9 0 20])
end


figure
nds='7654321'
clear t
for par=1:8
    subplot(2,4,par)
    for i=1:6
        tmp=param(par,i,:)<50;
        tmp=tmp(:);
        t(i,1)=nanmean(param(par,i,tmp));
        t(i,2)=nanstd(param(par,i,tmp))/sqrt(sum(~isnan(param(par,i,tmp))));
    end
    errorbar(t(:,1),t(:,2))
    title(m(par,:))    
    set(gca,'XTick',1:7,'xticklabel',{'ND7','ND6','ND5','ND4','ND3','ND2','ND1'})
%     axis([0 8 0 20])
end




% difference between ND2 and ND3 in amplitude vs in parameters
clear w
for j=1:6
    a=ampl(:,j+1)-ampl(:,j);
    for i=1:8
        tmp=reshape(param(i,j,:)<50&param(i,j+1,:)<50,1,53);
        b=reshape(param(i,j+1,tmp)-param(i,j,tmp),1,sum(tmp));
        if sum(tmp)>0
            w(j,i)=corr(a(tmp),b');
        end
    end
end
figure;imagesc(abs(w))
set(gca,'YTick',1:6,'yticklabel',{'ND7-ND6','ND6-5','ND5-4','ND4-3','ND3-2','ND2-1'})

figure
plot(abs(w'))
legend({'ND7-ND6','ND6-5','ND5-4','ND4-3','ND3-2','ND2-1'})


k=zeros(53,8);
for j=1:53
    for i=1:8        
        k(j,i)=corr(ampl(j,:)',param(i,:,j)');
    end
end
nanmean(k)







% relative to ND4
clear commRatio
for j=1:7
    for i=1:7
        switch i
            case 1
                a=spont(4,:)./spont(j,:);
            case 2
                a=LFC(4,:)./LFC(j,:);
            case 3
                a=MFC(4,:)./MFC(j,:);
            case 4
                a=HFC(4,:)./HFC(j,:);
            case 5
                a=LAC(4,:)./LAC(j,:);
            case 6
                a=MAC(4,:)./MAC(j,:);
            case 7
                a=HAC(4,:)./HAC(j,:);
            case 8
                a=tot(4,:)./tot(j,:);
        end
        a(isinf(a))=NaN;
        b(1:24,i,j)=a;
        commRatio(i,1,j)=nanmean(a);
        commRatio(i,2,j)=nanstd(a)/sqrt(sum(~isnan(a)));
    end
end
figure
nds='7654321'
for i=1:7
    subplot(2,4,i)
    errorbar(commRatio(:,1,i),commRatio(:,2,i))
    set(gca,'XTick',1:7,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    title(['ND4 vs ND', nds(i)])
end

figure
nds='7654321'
for j=1:7
    subplot(2,4,j)
    for i=1:7
        plot(ones(1,24)*i,b(:,i,j),'*')
        hold on
        plot(i,nanmean(b(:,i,j)),'r*')
    end
    line([0,8],[1,1],'color','k')
    axis([0,8,0,3])
    set(gca,'XTick',1:7,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    title(['ND4 vs ND', nds(j)])
end

% relative to ND1
clear commRatio
for j=1:7
    for i=1:7
        switch i
            case 1
                a=spont(7,:)./spont(j,:);
            case 2
                a=LFC(7,:)./LFC(j,:);
            case 3
                a=MFC(7,:)./MFC(j,:);
            case 4
                a=HFC(7,:)./HFC(j,:);
            case 5
                a=LAC(7,:)./LAC(j,:);
            case 6
                a=MAC(7,:)./MAC(j,:);
            case 7
                a=HAC(7,:)./HAC(j,:);
            case 8
                a=tot(7,:)./tot(j,:);
        end
        a(isinf(a))=NaN;
        b(1:24,i,j)=a;
        commRatio(i,1,j)=nanmean(a);
        commRatio(i,2,j)=nanstd(a)/sqrt(sum(~isnan(a)));
    end
end

figure
nds='7654321'
for j=1:7
    subplot(2,4,j)
    for i=1:7
        plot(ones(1,24)*i,b(:,i,j),'*')
        hold on
        plot(i,nanmean(b(:,i,j)),'r*')
    end
    line([0,8],[1,1],'color','k')
    axis([0,8,0,3])
    set(gca,'XTick',1:7,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    title(['ND1 vs ND', nds(j)])
end




figure
nds='7654321'
for j=1:7
    subplot(2,4,j)
    for i=1:24
        plot(1:7,b(i,:,j))
        hold on
%         plot(i,nanmean(b(:,i,j)),'r*')
    end
    line([0,8],[1,1],'color','k')
%     axis([0,8,0,3])
    set(gca,'XTick',1:7,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    title(['ND4 vs ND', nds(j)])
end




a=LAC(rang,:)./spont(rang,:)*100;
a(isinf(a))=100;
a(isnan(a))=0;
figure
errorbar(mean(a(:,3:end),2),std(a(:,3:end),0,2)/sqrt(24))

figure
scatterhist(a(4,:),a(6,:),'nbins',10)

for i=1:8
    switch i
        case 1
            a(1:7,1:24,i)=spont(rang,:)./spont(rang,:)*100;
        case 2
            a(1:7,1:24,i)=LFC(rang,:)./spont(rang,:)*100;
        case 3
            a(1:7,1:24,i)=MFC(rang,:)./spont(rang,:)*100;
        case 4
            a(1:7,1:24,i)=HFC(rang,:)./spont(rang,:)*100;
        case 5
            a(1:7,1:24,i)=LAC(rang,:)./spont(rang,:)*100;
        case 6
            a(1:7,1:24,i)=MAC(rang,:)./spont(rang,:)*100;
        case 7
            a(1:7,1:24,i)=HAC(rang,:)./spont(rang,:)*100;
        case 8
            a(1:7,1:24,i)=tot(rang,:)./spont(rang,:)*100;
    end
end
a(isinf(a))=100;
a(isnan(a))=0;

figure
nds='7654321'
for j=1:7
    subplot(2,4,j)
    for i=1:8
        plot(ones(1,24)*i,b(j,:,i),'*')
        hold on
        k=b(j,:,i)<5000;
        plot(i,mean(b(j,k,i)),'r*','MarkerSize',10)
    end
    line([0,8],[1,1],'color','k')
    axis([0,8,0,5000])
    set(gca,'XTick',1:7,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    title(['ND', nds(j)])
end



figure
nds='7654321'
for j=1:7
    subplot(2,4,j)
    for i=1:8
        plot(ones(1,24)*i,b(j,:,i),'*')
        hold on
        k=b(j,:,i)<5000;
        plot(i,mean(b(j,k,i)),'r*','MarkerSize',10)
    end
    line([0,8],[1,1],'color','k')
    axis([0,8,0,1500])
    set(gca,'XTick',1:7,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    title(['ND', nds(j)])
end



clear b
for j=1:7
    b(j,1:24,1:8)=a(4,:,:)./a(j,:,:);
end
b(isinf(b))=NaN;

figure
nds='7654321'
for j=1:7
    subplot(2,4,j)
    for i=1:8
        plot(ones(1,24)*i,b(j,:,i),'*')
        hold on
        plot(i,nanmean(b(j,:,i)),'r*')
    end
    line([0,8],[1,1],'color','k')
%     axis([0,8,0,3])
    set(gca,'XTick',1:8,'xticklabel',{'spo','LFC','MFC','HFC','LAC','MAC','HAC','tot'})
    title(['ND4 vs ND', nds(j)])
end



figure
nds='7654321'
for j=1:7
    subplot(2,4,j)
    for i=1:24
        plot(1:2,reshape(b(j,i,4:5),1,2))
        hold on
%         plot(i,nanmean(b(:,i,j)),'r*')
    end
    line([0,3],[1,1],'color','k')
%     axis([0,8,0,3])
    set(gca,'XTick',1:2,'xticklabel',{'HFC','LAC'})
    title(['ND4 vs ND', nds(j)])
end

for i=1:24
    tmp=max(abs(ampl(i,1:7)));
    c(i,1:7)=abs(ampl(i,1:7))/tmp;    
end

abs(ampl(:,[6 4 3]))
figure
for i=1:8
    subplot(2,4,i)
    plot(c,a(1:7,:,i)','.')
    legend(['7654321']')
    title([m(i,:)])
    axis([0,1,0,5000])
end
hold on
plot(a(1,:,4),ampl(:,1),'b.')

k=zeros(24,8);
for j=1:24
    for i=1:8
        k(j,i)=corr(c(j,1:7)',a(:,j,i));
    end
end
nanmean(k)