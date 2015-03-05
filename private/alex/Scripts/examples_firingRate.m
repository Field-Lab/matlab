cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121023'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end






figure
cnt=11;
st=96;
for i=1:2:5
    subplot(3,1,(i+1)/2)
    plot(mean(FiringRate(:,st+i:6:st+24,cnt),2))
    mean(mean(FiringRate(:,st+i:6:st+24,cnt),2))
    hold on
    plot(65001:130000,mean(FiringRate(:,st+i+1:6:st+24,cnt),2))
    mean(mean(FiringRate(:,st+i+1:6:st+24,cnt),2))
end


figure
hold on
cnt=11;
st=96;
a=correctedProtocols(1:362,2,st+1);
a(1)=30;
plot(correctedProtocols(1:362,1,st+1),a+180,'k')
plot(1:1984,ones(1984,1)*210,'k')
a=correctedProtocols(1:362,2,st+2);
a(1)=30;
plot(correctedProtocols(1:362,1,st+2)+8500,a+180,'k')
plot((1:1984)+8500,ones(1984,1)*210,'k')

plot(mean(FiringRate(1:8004,st+1:6:st+24,cnt),2),'b')
mean(mean(FiringRate(1:8004,st+1:6:st+24,cnt),2))

plot(8501:16504,mean(FiringRate(1:8004,st+2:6:st+24,cnt),2),'b')
mean(mean(FiringRate(1:8004,st+2:6:st+24,cnt),2))

title('High and Low contrast stimulus snippets and corresponding firing rate of a ganglion cell')
text(3000,240,'high contrast stimulus, a.u.')
text(11000,240,'low contrast stimulus,a.u.')
text(3000,150,'firing rate, Hz')
text(11000,150,'firing rate, Hz')
xlabel('time, ms')



figure
hold off
cnt=13;
st=96;
a=correctedProtocols(1:362,2,st+1);
a(1)=30;
plot(correctedProtocols(1:362,1,st+1),a/6+25,'k')
hold on
plot(1:1984,ones(1984,1)*30,'k')
a=correctedProtocols(1:362,2,st+2);
a(1)=30;
plot(correctedProtocols(1:362,1,st+2)+8500,a/6+25,'k')
plot((1:1984)+8500,ones(1984,1)*30,'k')

plot(mean(FiringRate(1:8004,st+1:6:st+24,cnt),2),'b')
mean(mean(FiringRate(1:8004,st+1:6:st+24,cnt),2))

plot(8501:16504,mean(FiringRate(1:8004,st+2:6:st+24,cnt),2),'b')
mean(mean(FiringRate(1:8004,st+2:6:st+24,cnt),2))

title('High and Low contrast stimulus snippets and corresponding firing rate of a ganglion cell')
text(3000,34,'high contrast stimulus, a.u.')
text(11000,34,'low contrast stimulus,a.u.')
text(3000,19,'firing rate, Hz')
text(11000,19,'firing rate, Hz')
xlabel('time, ms')




figure
hold off
cnt=13;
st=96;
a=correctedProtocols(1:3000,2,st+1);
a(1)=30;
plot(correctedProtocols(1:3000,1,st+1),a/6+25,'color',[0.7 0.7 0.7])
hold on
plot(1:1984,ones(1984,1)*30,'color',[0.7 0.7 0.7])
a=correctedProtocols(1:3000,2,st+2);
a(1)=30;
plot(correctedProtocols(1:3000,1,st+2)+57000,a/6+25,'color',[0.7 0.7 0.7])
plot((1:1984)+57000,ones(1984,1)*30,'color',[0.7 0.7 0.7])

plot(mean(FiringRate(1:52000,st+1:6:st+24,cnt),2),'k')

m0=mean(mean(FiringRate(1:2000,st+1:6:st+24,cnt),2))
va0=std(mean(FiringRate(1:2000,st+1:6:st+24,cnt),2))

m=mean(mean(FiringRate(2000:52000,st+1:6:st+24,cnt),2))
va=std(mean(FiringRate(2000:52000,st+1:6:st+24,cnt),2))

plot(57001:109000,mean(FiringRate(1:52000,st+2:6:st+24,cnt),2),'k')
mean(mean(FiringRate(1:52000,st+2:6:st+24,cnt),2))
m1=mean(mean(FiringRate(2000:52000,st+2:6:st+24,cnt),2))
va1=std(mean(FiringRate(2000:52000,st+2:6:st+24,cnt),2))


tmp=mean(LinearFilter(:,1,st+1:6:st+24,cnt),3);
tmp=tmp/sum(abs(tmp));
plot(1:104:52000,tmp*500-12,'k','linewidth',2)
hold on
tmp=mean(LinearFilter(:,1,st+2:6:st+24,cnt),3);
tmp=tmp/sum(abs(tmp));
plot((1:104:52000)+57000,tmp*500-12,'k','linewidth',2)

bar(53500,m0,800,'facecolor','r')
bar(54500,m,800,'facecolor',[0.3 0.8 0.1])
bar(55500,m1,800,'facecolor','b')
errorbar(53500:1000:55500,[m0 m m1],[va0 va va1],'.k','linewidth',2)
text(48000,-1,'spontaneous firing rate mean and std','color','r')
text(48000,-2,'high contrast firing rate mean and std','color',[0.3 0.8 0.1])
text(48000,-3,'low contrast firing rate mean and std','color','b')

title('High and Low contrast stimulus snippets (50s) and corresponding firing rate of a ganglion cell')
text(3000,37,'high contrast stimulus, a.u.')
text(59000,37,'low contrast stimulus,a.u.')
text(3000,17,'firing rate, Hz')
text(59000,22,'firing rate, Hz')
text(12000,-8,'sta,a.u.')
text(70000,-8,'sta,a.u.')
axis([0 110000 -20 40])