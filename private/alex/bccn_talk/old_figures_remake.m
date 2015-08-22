load('/Users/alexth/Desktop/old_stuff/my_data/analysis/summary_all.mat')


figure
nds='87654321'
data=zeros(4500,8,2);
for i=1:8
    subplot(2,4,i)
    tmp=wf{i}(:,129);
    hold off
    plot(tmp,'color','k','linewidth',2)
    data(:,i,1)=tmp;
    
    hold on
    tmp=bf{i}(:,129);
    plot(4701:9200,tmp,'color','k','linewidth',2)
    data(:,i,2)=tmp;
    
    axis([1 9200 0 100])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55]+30,'color','k')
    line([500 2500],[60,60]+30,'color','k')
    line([2500 4500],[55,55]+30,'color','k')
    line([500 500],[55,60]+30,'color','k')
    line([2500 2500],[55,60]+30,'color','k')
    
    line([1 500]+4700,[55,55]+30,'color','k')
    line([500 2500]+4700,[50,50]+30,'color','k')
    line([2500 4500]+4700,[55,55]+30,'color','k')
    line([500 500]+4700,[50,55]+30,'color','k')
    line([2500 2500]+4700,[50,55]+30,'color','k')


    title(['ND',nds(i)])
end


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=201;
beg2=4901;
add=10;
ndlabel=10;
stim_height=10;

for i=1:8
    
    minn=abs(min([data(:,i,1); data(:,i,2)]))+add;
    maxx=max([data(:,i,1); data(:,i,2)]);
    zeroCoord(i)=minn+const;
    if i==1
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,i,1)+minn+const,'k','linewidth',2);    
    plot(beg2:beg2+4499,data(:,i,2)+minn+const,'k','linewidth',2);
    
%     text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg2+4500 -35 const])

k=-20;
for i=1:4
    set(hline(i),'ydata',[k zeroCoord(8)+maxx])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500],[stim_height,stim_height]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500],[-1,stim_height+1]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500],[-1,stim_height+1]+k,'color','k','linewidth',3);

line([beg1 beg1+500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500]+beg2-beg1,[-stim_height,-stim_height]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500]+beg2-beg1,[-1.5-stim_height,1.5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500]+beg2-beg1,[-1.5-stim_height,1.5]+k,'color','k','linewidth',3);

line([beg2-beg1-500 beg2-beg1-500],[zeroCoord(2)+20,zeroCoord(2)+40],'color','k','linewidth',3);
% text(beg2-beg1-500+100,zeroCoord(2)+30,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])


