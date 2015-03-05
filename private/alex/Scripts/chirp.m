clear
cd('S:\user\alexandra\scripts')
date='20120928';
path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

% path2save=[path2data,'Full_info_5chirp\'];
% if ~exist(path2save,'file')
%     mkdir(path2save);
% end

chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);
flips=spike_info.flip_times{1}(:,1); % flips in ms, rounded
pxls=spike_info.flip_times{1}(:,2); % pxs in ms, normalized -1:1
pxls(pxls==-1)=30;
tmp1=zeros(flips(end)+1,1); % 1 column - color, 2 column - spikes
tmp1([1; flips(1:end-1)],1)=pxls;

for i=1:50
    subscr=(flips(1:end-1)+i);
    subscr(tmp1(flips(1:end-1)+i,1)~=0)=[];
    tmp1(subscr,1)=tmp1(subscr-1,1);
end
tmp1(tmp1==0)=30;
noStim=200:2500;
lowFreqChirp=3030:4935;
% midFreqChirp=4935:7953;
% midFreqChirp=4935:9720;
midFreqChirp=4935:10000;
highFreqChirp=10000:10970;
lowAmplChirp=11670:12600;
midAmplChirp=12600:16100;
highAmplChirp=16100:19600;
% figure
clear ampl
for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    cnt=1;
    for i=10:12:112
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
    for i=1:5:trialsChirp
        tmp(cnt,1:20000)=mean(chirp(i:i+4,1:20000));
        cnt=cnt+1;
    end
%     subplot(5,6,unit)
%     plot(tmp(5,:))
    for j=1:trialsChirp/5
        baseLine=mean(tmp(j,200:2200));
        based=abs(tmp(j,:)-baseLine);
        spont(j,unit)=sum(based(noStim))/length(noStim);
        LFC(j,unit)=sum(based(lowFreqChirp))/length(lowFreqChirp);
        MFC(j,unit)=sum(based(midFreqChirp))/length(midFreqChirp);
        HFC(j,unit)=sum(based(highFreqChirp))/length(highFreqChirp);
        LAC(j,unit)=sum(based(lowAmplChirp))/length(lowAmplChirp);
        MAC(j,unit)=sum(based(midAmplChirp))/length(midAmplChirp);
        HAC(j,unit)=sum(based(highAmplChirp))/length(highAmplChirp);
        tot(j,unit)=sum(tmp(j,:))/length(tmp);
    end
    % end
    
    % plot(tmp')
    
    % contrastMod=11650:19590;
    %
    % for trial=1:length(tmp)
    %     a=tmp(1,contrastMod)-mean(tmp(1,1:1500));
    %
    % end
    %  b=tmp(1,contrastMod)-mean(tmp(1,1:1500));
    %  plot(b)
    %  plot(chirp(1:5,contrastMod)')
    %  hold on
    %   plot(b,'linewidth',3)
    %   plot(tmp1(contrastMod)+20,'r','linewidth',3)
    %
%     
%     [~,ind]=findpeaks(tmp1(contrastMod));
%     ind2=ind+250;
%     clear k l
%     for j=1:14
%         b=tmp(j,contrastMod)-mean(tmp(j,1:1500));
%         c=[zeros(1,38) b(ind(1)-250:ind(1)+750)>0 zeros(1,6902)];
%         c=logical(c);
%         d=[zeros(1,6574) b(ind(end-2)-250:ind(end-2)+750)>0 zeros(1,366)];
%         d=logical(d);
%         if isempty(c)
%             k(j,unit)=NaN;
%         else
%             k(j,unit)=sum(b(c))/sum(b(d));
%         end
%   
%     end
    % nds='8765432123456'
    % clear mins maxs
%     for j=1:12
%         %     subplot(3,4,j)
%         %     bar(k(:,:,j))
%         mins(j,1:4,unit)=k(1:4,1,j)/mean(k(end-3:end,1,j));
%         maxs(j,1:4,unit)=k(1:4,2,j)/mean(k(end-3:end,2,j));
%         %     title(['ND',nds(j),'  ','mean ',int2str(round(l(j)))])
%         %     axis([0 15 -50 100])
%     end
end
% subplot(5,6,30)
% plot(tot(5,:))


figure
rang=1:7;
commMean(1:7,1)=mean(spont(rang,:),2);
commMean(1:7,2)=mean(LFC(rang,:),2);
commMean(1:7,3)=mean(MFC(rang,:),2);
commMean(1:7,4)=mean(HFC(rang,:),2);
commMean(1:7,5)=mean(LAC(rang,:),2);
commMean(1:7,6)=mean(MAC(rang,:),2);
commMean(1:7,7)=mean(HAC(rang,:),2);
commMean(1:7,8)=mean(tot(rang,:),2);

commStd(1:7,1)=std(spont(rang,:),0,2)/sqrt(size(spont,2));
commStd(1:7,2)=std(LFC(rang,:),0,2)/sqrt(size(spont,2));
commStd(1:7,3)=std(MFC(rang,:),0,2)/sqrt(size(spont,2));
commStd(1:7,4)=std(HFC(rang,:),0,2)/sqrt(size(spont,2));
commStd(1:7,5)=std(LAC(rang,:),0,2)/sqrt(size(spont,2));
commStd(1:7,6)=std(MAC(rang,:),0,2)/sqrt(size(spont,2));
commStd(1:7,7)=std(HAC(rang,:),0,2)/sqrt(size(spont,2));
commStd(1:7,8)=std(tot(rang,:),0,2)/sqrt(size(spont,2));

 m=['spo';'LFC';'MFC';'HFC';'LAC';'MAC';'HAC';'tot']
for i=1:8
    subplot(2,4,i)
    hold on
    errorbar(commMean(:,i),commStd(:,i),'b')
    set(gca,'XTick',1:7,'xticklabel',{'ND7','ND6','ND5','ND4','ND3','ND2','ND1'})
    title([m(i,:)])%,'   ALL red 1002 blue 0928'])
end

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
rectangle('position',[0,0,22500,60],'facecolor',[0.8 0.8 0.8])
rectangle('position',[noStim(1),0,noStim(end)-noStim(1),60],'facecolor',[0.8 0.6 0.6])
rectangle('position',[lowFreqChirp(1),0,lowFreqChirp(end)-lowFreqChirp(1),60],'facecolor',[0.7 0.9 0.7])
rectangle('position',[midFreqChirp(1),0,midFreqChirp(end)-midFreqChirp(1),60],'facecolor',[0.6 0.8 0.6])
rectangle('position',[highFreqChirp(1),0,highFreqChirp(end)-highFreqChirp(1),60],'facecolor',[0.5 0.7 0.5])
rectangle('position',[lowAmplChirp(1),0,lowAmplChirp(end)-lowAmplChirp(1),60],'facecolor',[0.7 0.9 0.7])
rectangle('position',[midAmplChirp(1),0,midAmplChirp(end)-midAmplChirp(1),60],'facecolor',[0.6 0.8 0.6])
rectangle('position',[highAmplChirp(1),0,highAmplChirp(end)-highAmplChirp(1),60],'facecolor',[0.5 0.7 0.5])
hold on
plot(tmp1,'k')
axis([0,22500,0,60])


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