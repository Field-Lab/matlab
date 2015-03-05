%% Fig.E
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_E')

figure
set(gcf,'position',[1784         240         446         715])

p=subplot('position',[0.01 0.01 0.9 0.92]);
hold on
const=0;
beg1=201;
beg2=4901;
add=2;

for i=1:8
    
    minn=abs(min([data(:,1,i); data(:,2,i); data(:,3,i); data(:,4,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i); data(:,3,i); data(:,4,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        plot(beg1:beg1+4499,data(:,1,i)+minn+const,'r','linewidth',2);
        plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
        mleg=legend({'non-typical ON cells','typical ON cells'});
        legend('boxoff')
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,1,i)+minn+const,'r','linewidth',2);
    plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
    
    plot(beg2:beg2+4499,data(:,3,i)+minn+const,'r','linewidth',2);
    plot(beg2:beg2+4499,data(:,4,i)+minn+const,'b','linewidth',2);
    
    text(2,zeroCoord(i)+10,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end
set(mleg,opt_title,'position',[0.2 0.91 0.7 0.05])
set(mleg,'fontsize',10)
const=const+15;
axis([1 beg2+4500 -20 const])

k=-10;
for i=1:4
    set(hline(i),'ydata',[k zeroCoord(8)+maxx])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500],[5,5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500],[-1,6]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500],[-1,6]+k,'color','k','linewidth',3);

line([beg1 beg1+500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500]+beg2-beg1,[-5,-5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500]+beg2-beg1,[-6.5,1.5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500]+beg2-beg1,[-6.5,1.5]+k,'color','k','linewidth',3);




line([beg2-beg1-500 beg2-beg1-500],[zeroCoord(2)+10,zeroCoord(2)+30],'color','k','linewidth',3);
text(beg2-beg1-500+100,zeroCoord(2)+20,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\E.emf')



%% Fig.X
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_X.mat')


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=201;
beg2=4901;
add=70;
ndlabel=60;
stim_height=15;

for i=1:8
    
    minn=abs(min([data(:,1,i); data(:,2,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,2,i)+minn+const,'k','linewidth',2);    
    plot(beg2:beg2+4499,data(:,1,i)+minn+const,'k','linewidth',2);

    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg2+4500 -10 const])

k=10;
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

line([beg2-beg1-500 beg2-beg1-500],[zeroCoord(2)+60,zeroCoord(2)+120],'color','k','linewidth',3);
text(beg2-beg1-500+100,zeroCoord(2)+90,'60Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\X.emf')



%% Fig.Y
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Y.mat')


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=201;
beg2=4901;
add=70;
ndlabel=30;
stim_height=15;

for i=1:8
    
    minn=abs(min([data(:,1,i); data(:,2,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,2,i)+minn+const,'k','linewidth',2);    
    plot(beg2:beg2+4499,data(:,1,i)+minn+const,'k','linewidth',2);

    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg2+4500 -30 const])

k=-10;
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




line([beg2-beg1-500 beg2-beg1-500],[zeroCoord(2)+40,zeroCoord(2)+80],'color','k','linewidth',3);
text(beg2-beg1-500+100,zeroCoord(2)+60,'40Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])


saveas(gcf,'S:\user\alexandra\thesis\figures\Y.emf')



%% Fig.Z
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z.mat')


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=201;
beg2=4901;
add=10;
ndlabel=10;
stim_height=5;

for i=1:8
    
    minn=abs(min([data(:,1,i); data(:,2,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,1,i)+minn+const,'k','linewidth',2);    
    plot(beg2:beg2+4499,data(:,2,i)+minn+const,'k','linewidth',2);

    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg2+4500 -15 const])

k=-5;
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

line([beg2-beg1-500 beg2-beg1-500],[zeroCoord(2)+10,zeroCoord(2)+30],'color','k','linewidth',3);
text(beg2-beg1-500+100,zeroCoord(2)+20,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\Z.emf')



%% Fig.Z9
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z9.mat')


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
    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
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
text(beg2-beg1-500+100,zeroCoord(2)+30,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\Z9.emf')


%% Fig.Z10
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z10.mat')


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=201;
beg2=4901;
add=10;
ndlabel=10;
stim_height=5;

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
    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg2+4500 -20 const])

k=-10;
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
text(beg2-beg1-500+100,zeroCoord(2)+30,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])


saveas(gcf,'S:\user\alexandra\thesis\figures\Z10.emf')


%% Fig.Z15
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z15.mat')


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=201;
beg2=4901;
add=5;
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
    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg2+4500 -40 const])

k=-25;
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
text(beg2-beg1-500+100,zeroCoord(2)+30,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])


saveas(gcf,'S:\user\alexandra\thesis\figures\Z15.emf')


%% Fig.R
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R.mat')


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=201;
beg2=4901;
add=5;
ndlabel=5;
stim_height=5;

for i=1:8
    
    minn=abs(min([data(:,1,i); data(:,2,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,1,i)+minn+const,'k','linewidth',2);    
    plot(beg2:beg2+4499,data(:,2,i)+minn+const,'k','linewidth',2);

    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg2+4500 -20 const])

k=-10;
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

line([beg2-beg1-500 beg2-beg1-500],[zeroCoord(2)+10,zeroCoord(2)+30],'color','k','linewidth',3);
text(beg2-beg1-500+100,zeroCoord(2)+20,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\R.emf')



%% Fig.R3
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R3_R4_R5_R6.mat','data','bardata')


figure
set(gcf,'position',[1784         240         446         715])
p=subplot('position',[0.01 0.01 0.9 0.97]);
hold on
const=0;
beg1=101;
add=5;
ndlabel=10;
stim_height=5;

for i=1:8

    minn=abs(min([data(:,1,i); data(:,2,i); data(:,3,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i); data(:,3,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        plot(beg1:beg1+4499,data(:,1,i)+minn+const,'color',[0.3 0.8 0.1],'linewidth',2);
        plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
        plot(beg1:beg1+4499,data(:,3,i)+minn+const,'r','linewidth',2);
        mleg=legend({'transient','gap-like','sustained'});
        legend('boxoff')
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end
    
    line([0 beg1+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,1,i)+minn+const,'color',[0.3 0.8 0.1],'linewidth',2);
    plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
    plot(beg1:beg1+4499,data(:,3,i)+minn+const,'r','linewidth',2);
   
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg1+4500 -20 const])

k=-10;
for i=1:2
    set(hline(i),'ydata',[k zeroCoord(8)+maxx])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500],[stim_height,stim_height]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500],[-1,stim_height+1]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500],[-1,stim_height+1]+k,'color','k','linewidth',3);


line([beg1+3800 beg1+3800],[zeroCoord(2)+10,zeroCoord(2)+30],'color','k','linewidth',3);
text(beg1+3900,zeroCoord(2)+20,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])


saveas(gcf,'S:\user\alexandra\thesis\figures\R3.emf')


%% Fig.P
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_P.mat')
figure
set(gcf,'position',[1784         240         446         715])

p=subplot('position',[0.01 0.01 0.9 0.95]);
hold on
const=0;
beg1=201;
beg2=4901;
add=5;
ndlabel=10;
stim_height=5;

for i=1:8
    
    minn=abs(min([data(:,1,i); data(:,2,i); data(:,3,i); data(:,4,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i); data(:,3,i); data(:,4,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        plot(beg1:beg1+4499,data(:,1,i)+minn+const,'r','linewidth',2);
        plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
        mleg=legend({'non-typical ON cells','typical ON cells'});
        legend('boxoff')
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,1,i)+minn+const,'r','linewidth',2);
    plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
    
    plot(beg2:beg2+4499,data(:,3,i)+minn+const,'r','linewidth',2);
    plot(beg2:beg2+4499,data(:,4,i)+minn+const,'b','linewidth',2);
    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end
set(mleg,opt_title,'position',[0.2 0.94 0.7 0.05])
set(mleg,'fontsize',10)
const=const+15;
axis([1 beg2+4500 -20 const])

k=-10;
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


line([beg2-beg1-500 beg2-beg1-500],[zeroCoord(2)+10,zeroCoord(2)+30],'color','k','linewidth',3);
text(beg2-beg1-500+100,zeroCoord(2)+20,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\P.emf')


%% Fig.R5
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R3_R4_R5_R6.mat','data','bardata')

data=data(:,:,3:5);
figure
set(gcf,'position',[2112         404         302         480])
p=subplot('position',[0.01 0.01 0.93 0.99]);
hold on
const=0;
beg1=101;
add=5;
ndlabel=10;
stim_height=5;

for i=1:3

    minn=abs(min([data(:,i,1); data(:,i,2); data(:,i,3)]))+add;
    maxx=max([data(:,i,1); data(:,i,2); data(:,i,3)]);
    zeroCoord(i)=minn+const;
    if i==1
        plot(beg1:beg1+4499,data(:,i,1)+minn+const,'color',[0.3 0.8 0.1],'linewidth',2);
        plot(beg1:beg1+4499,data(:,i,2)+minn+const,'b','linewidth',2);
        plot(beg1:beg1+4499,data(:,i,3)+minn+const,'r','linewidth',2);
        mleg=legend({'ND6','ND5','ND4'});
        legend('boxoff')
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end
    
    line([0 beg1+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,i,1)+minn+const,'color',[0.3 0.8 0.1],'linewidth',2);
    plot(beg1:beg1+4499,data(:,i,2)+minn+const,'b','linewidth',2);
    plot(beg1:beg1+4499,data(:,i,3)+minn+const,'r','linewidth',2);
    
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

text(3000,zeroCoord(1)+20,'transient',opt_text,'fontsize',14)
text(3400,zeroCoord(2)+20,'gap',opt_text,'fontsize',14)
text(3000,zeroCoord(3)+20,'sustained',opt_text,'fontsize',14)

const=const+15;
axis([1 beg1+4500 -10 const])

k=-5;
for i=1:2
    set(hline(i),'ydata',[k zeroCoord(3)+maxx])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500],[stim_height,stim_height]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500],[-1,stim_height+1]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500],[-1,stim_height+1]+k,'color','k','linewidth',3);


line([beg1+1500 beg1+1500],[zeroCoord(1)+10,zeroCoord(1)+30],'color','k','linewidth',3);
text(beg1+1600,zeroCoord(1)+20,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\R5.emf')


%% Fig.R11
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R11.mat','data','bardata')


figure
set(gcf,'position',[1784         240         297         715])
p=subplot('position',[0.01 0.01 0.98 0.97]);
hold on
const=0;
beg1=1;
add=5;
ndlabel=10;
stim_height=7;

for i=1:8

    minn=abs(min([data(:,1,i); data(:,2,i); data(:,3,i)]))+add;
    maxx=max([data(:,1,i); data(:,2,i); data(:,3,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        plot(beg1:beg1+4499,data(:,1,i)+minn+const,'color',[0.3 0.8 0.1],'linewidth',2);
        plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
        plot(beg1:beg1+4499,data(:,3,i)+minn+const,'r','linewidth',2);
        mleg=legend({'transient','gap-like','sustained'});
        set(mleg,'position',[0.3 0.9 0.1 0.1])
        legend('boxoff')
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end
    
    line([0 beg1+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,data(:,1,i)+minn+const,'color',[0.3 0.8 0.1],'linewidth',2);
    plot(beg1:beg1+4499,data(:,2,i)+minn+const,'b','linewidth',2);
    plot(beg1:beg1+4499,data(:,3,i)+minn+const,'r','linewidth',2);
   
    text(2,zeroCoord(i)+ndlabel,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
      
end

const=const+15;
axis([1 beg1+4500 -20 const])

k=-10;
for i=1:2
    set(hline(i),'ydata',[k zeroCoord(8)+maxx])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500],[-stim_height,-stim_height]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500],[-1.5-stim_height,1.5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500],[-1.5-stim_height,1.5]+k,'color','k','linewidth',3);


line([beg1+1800 beg1+1800],[zeroCoord(2)+10,zeroCoord(2)+30],'color','k','linewidth',3);
text(beg1+1900,zeroCoord(2)+20,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

saveas(gcf,'S:\user\alexandra\thesis\figures\R11.emf')



%% Fig.Z3
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z3.mat')

figure
set(gcf,'position',[2126         443         651         382])
nds='666555444333222';

corn_posy=0.03;
corn_posx=0.07:0.155:1;
size_x=0.14;
size_y=0.96;
set(gcf,'color',[1 1 1])
beg1=1;
add=2;
ndlabel=10;
stim_height=7;
k=-10;
for j=1:6
    subplot('position',[corn_posx(j),corn_posy,size_x,size_y])

    hold on
    const=0;
    
    for i=1:3
        minn=abs(min(min(data(:,i,:))))+add;
        maxx=max(max(data(:,i,:)));
        zeroCoord(i)=minn+const;
        if i==1
            hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
            hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        end
        
        line([0 beg1+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
        plot(data(:,i,j)+minn+const,'k','linewidth',2)
        
        if j==1
            text(-1500,zeroCoord(i)+ndlabel,['ND',int2str(7-i)],opt_text)
        end
        
        const=const+minn+maxx;
    end
     
    const=const+5;
    axis([1 beg1+4500 -20 const])

   
    for i=1:2
        set(hline(i),'ydata',[k zeroCoord(3)+maxx])
    end
    
    line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
    line([beg1+500 beg1+2500],[-stim_height,-stim_height]+k,'color','k','linewidth',3);
    line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
    line([beg1+500 beg1+500],[-1.5-stim_height,1.5]+k,'color','k','linewidth',3);
    line([beg1+2500 beg1+2500],[-1.5-stim_height,1.5]+k,'color','k','linewidth',3);
    set(gca,'xtick',0,'xticklabel','','ytick',0,'yticklabel','')
    
    if j==1
        line([beg1+2900 beg1+2900],[zeroCoord(1)+10,zeroCoord(1)+30],'color','k','linewidth',3);
        text(beg1+3100,zeroCoord(1)+20,'20Hz',opt_text)
    end
    set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])

end

    
saveas(gcf,'S:\user\alexandra\thesis\figures\Z3.emf')



%% Fig.Q
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Q')

figure
set(gcf,'position',[2053         151         658         715])

p=subplot('position',[0.01 0.01 0.45 0.92]);
hold on
const=0;
beg1=201;
beg2=4901;
add=2;

for i=1:8
    
    minn=abs(min([dataOfffirst(:,1,i); dataOfffirst(:,2,i);dataOfflast(:,1,i);dataOfflast(:,2,i)]))+add;
    maxx=max([dataOfffirst(:,1,i); dataOfffirst(:,2,i);dataOfflast(:,1,i);dataOfflast(:,2,i)]);
    zeroCoord(i)=minn+const;
    if i==1
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,dataOfffirst(:,1,i)+minn+const,'r','linewidth',2);
    plot(beg1:beg1+4499,dataOfflast(:,1,i)+minn+const,'b','linewidth',2);
    
    plot(beg2:beg2+4499,dataOfffirst(:,2,i)+minn+const,'r','linewidth',2);
    plot(beg2:beg2+4499,dataOfflast(:,2,i)+minn+const,'b','linewidth',2);
    
%     text(2,zeroCoord(i)+10,['ND',int2str(9-i)],opt_text)
    const=const+minn+maxx;
    const1(i)=const;
      
end

const=const+15;
axis([1 beg2+4500 -20 const])

k=-10;
for i=1:4
    set(hline(i),'ydata',[k zeroCoord(8)+maxx])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500],[5,5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500],[-1,6]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500],[-1,6]+k,'color','k','linewidth',3);

line([beg1 beg1+500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500]+beg2-beg1,[-5,-5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500]+beg2-beg1,[-6.5,1.5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500]+beg2-beg1,[-6.5,1.5]+k,'color','k','linewidth',3);




line([beg2+3000 beg2+3000],[zeroCoord(2)+20,zeroCoord(2)+40],'color','k','linewidth',3);
text(beg2+3100,zeroCoord(2)+30,'20Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])




p=subplot('position',[0.53 0.01 0.45 0.92]);
hold on
const=0;
beg1=201;
beg2=4901;
add=2;

for i=1:8
    
    minn=abs(min([dataOnfirst(:,1,i); dataOnfirst(:,2,i);dataOnlast(:,1,i);dataOnlast(:,2,i)]))+add;
    maxx=max([dataOnfirst(:,1,i); dataOnfirst(:,2,i);dataOnlast(:,1,i);dataOnlast(:,2,i)]);
%     zeroCoord(i)=minn+const;
    if i==1
        plot(beg1:beg1+4499,dataOnfirst(:,1,i)+zeroCoord(i),'r','linewidth',2);
        plot(beg1:beg1+4499,dataOnlast(:,1,i)+zeroCoord(i),'b','linewidth',2);
        mleg=legend({'first trial','last trial'});
        legend('boxoff')
        hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
        hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
    end    
    
    line([0 beg2+4500],[zeroCoord(i) zeroCoord(i)],'color',[1 1 1]*0.5,'linewidth',2);
    
    plot(beg1:beg1+4499,dataOnfirst(:,1,i)+zeroCoord(i),'r','linewidth',2);
    plot(beg1:beg1+4499,dataOnlast(:,1,i)+zeroCoord(i),'b','linewidth',2);
    
    plot(beg2:beg2+4499,dataOnfirst(:,2,i)+zeroCoord(i),'r','linewidth',2);
    plot(beg2:beg2+4499,dataOnlast(:,2,i)+zeroCoord(i),'b','linewidth',2);
    
    text(-1000,zeroCoord(i)+10,['ND',int2str(9-i)],opt_text)
%     const=const+minn+maxx;
      const=const1(i);
end
set(mleg,opt_title,'position',[0.4 0.93 0.2 0.05])
set(mleg,'fontsize',10)
const=const+15;
axis([1 beg2+4500 -20 const])

k=-10;
for i=1:4
    set(hline(i),'ydata',[k zeroCoord(8)+maxx])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500],[5,5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500],[-1,6]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500],[-1,6]+k,'color','k','linewidth',3);

line([beg1 beg1+500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+2500]+beg2-beg1,[-5,-5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+4500]+beg2-beg1,[0,0]+k,'color','k','linewidth',3);
line([beg1+500 beg1+500]+beg2-beg1,[-6.5,1.5]+k,'color','k','linewidth',3);
line([beg1+2500 beg1+2500]+beg2-beg1,[-6.5,1.5]+k,'color','k','linewidth',3);

set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])


saveas(gcf,'S:\user\alexandra\thesis\figures\Q.emf')