a=dir('S:\data\alexandra\MEA_data\20121216\Filters\*HC*.mat')
HighFilters=zeros(500,size(a,1));
for j=1:size(a,1)
    load(['S:\data\alexandra\MEA_data\20121216\Filters\',a(j).name])
    HighFilters(1:500,j)=mean(filters(:,[3,4,8,13,18,27]),2)/10e3;
end

figure
cnt=1;

nds='87654321234';
for i=1:12:size(a,1)    
    subplot(3,4,cnt)
    line([0 500],[0,0],'color','k')
    hold on
    goto=min(i+11,size(a,1));
    plot(HighFilters(1:500,i:3:goto),'b')
    plot(HighFilters(1:500,i+1:3:goto),'r')
    plot(HighFilters(1:500,i+2:3:goto),'k')
%     axis([0, 500, -160, 70])
    title(['ND',nds(cnt)])
    set(gca,'xtick',0,'xticklabel','')%,'ytick',0,'xticklabel','')%[round(min(mean(filters,2)/10e3)*100)/100,0,round(max(mean(filters,2)/10e3)*100)/100])
    cnt=cnt+1;
end


a=dir('S:\data\alexandra\MEA_data\20121216\Filters\*LCseq*.mat')
LowFilters=zeros(500,132);
for j=1:132
    load(['S:\data\alexandra\MEA_data\20121216\Filters\',a(j).name])
    LowFilters(1:500,j)=mean(filters(:,31:60),2)/10e3;
end


cnt=1;

nds='87654321234';
for i=1:12:132    
    subplot(3,4,cnt)
    line([0 500],[0,0],'color','k')
    hold on
    plot(LowFilters(1:500,i:3:i+11),'b')
    plot(LowFilters(1:500,i+1:3:i+11),'r')
    plot(LowFilters(1:500,i+2:3:i+11),'k')
    axis([0, 500, -10, 4])
    title(['ND',nds(cnt)])
    set(gca,'xtick',0,'xticklabel','')%,'ytick',0,'xticklabel','')%[round(min(mean(filters,2)/10e3)*100)/100,0,round(max(mean(filters,2)/10e3)*100)/100])
    cnt=cnt+1;
end