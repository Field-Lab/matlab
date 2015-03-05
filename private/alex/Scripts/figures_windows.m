cd('S:\user\alexandra\scripts')
clear
opt_axis=struct('fontsize',10,'fontname','arial');
opt_title=struct('fontweight','bold','fontsize',12,'fontname','arial');
opt_leg=struct('fontsize',12,'fontname','arial');
opt_label=struct('fontsize',10,'fontweight','bold','fontname','arial');
opt_text=struct('fontsize',10,'fontname','arial','fontweight','bold');
plnmr=[4 7 2 5 8 3 6 9];
corn_posy=[0.395 0.08 0.71 0.395 0.08 0.71 0.395 0.08]+0.01;
corn_posx=[0.08 0.08 0.405 0.405 0.405 0.73 0.73 0.73]+0.01;
size_x=0.25;
size_y=0.23;
sizes_1=[0.0600    0.7200    0.2870    0.2500];


save('S:\data\alexandra\MEA_data\analysis\plot_opts','opt_axis','opt_title','opt_leg','opt_label','opt_text','plnmr',...
    'corn_posy','corn_posx','size_x','size_y','sizes_1');

%% Figure M1
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\summaryzc_late','zc_HC','peak_HC','onOff')

figure
set(gcf,'position',[1962         503         600         404])

for i=1:8
    a(i)=mean(peak_HC(i,peak_HC(i,:)>0&onOff>0));
    b(i)=std(peak_HC(i,peak_HC(i,:)>0))/sqrt(sum(peak_HC(i,:)>0&onOff>0));
    c(i)=mean(peak_HC(i,peak_HC(i,:)>0&onOff<0));
    d(i)=std(peak_HC(i,peak_HC(i,:)>0))/sqrt(sum(peak_HC(i,:)>0&onOff<0));
end
errorbar(a,b,'r.')
hold on
errorbar(c,d,'b.')
clear a
for i=1:8
    a(i)=mean(zc_HC(i,zc_HC(i,:)>0&onOff>0));
    b(i)=std(zc_HC(i,zc_HC(i,:)>0))/sqrt(sum(zc_HC(i,:)>0&onOff>0));
    c(i)=mean(zc_HC(i,zc_HC(i,:)>0&onOff<0));
    d(i)=std(zc_HC(i,zc_HC(i,:)>0))/sqrt(sum(zc_HC(i,:)>0&onOff<0));
end 
errorbar(a,b,'m.')
hold on
errorbar(c,d,'c.')

m=legend({'peak ON cells','peak OFF cells','zero crossing ON cells','zero crossing OFF cells'});
set(m,opt_leg);
axis([0 9 70 220])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Latency, ms');
set(m,opt_label);
saveas(gcf,'S:\user\alexandra\thesis\figures\M1.emf')

%% Fig.M2 temporal development of latency
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_M2')

figure
set(gcf,'position',[1833         425        1119         505])
errorbar(coord,peak_on,peak_on_std,'r.')
hold on
errorbar(coord,peak_off,peak_off_std,'b.')
errorbar(coord,zc_on,zc_on_std,'m.')
errorbar(coord,zc_off,zc_off_std,'c.')

m=legend({'peak ON','peak OFF','zero crossing ON','zero crossing OFF'});
set(m,opt_leg);
axis([0 160 70 230])
set(gca,'xtick',9.5:20:160,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Latency, ms');
set(m,opt_label);
saveas(gcf,'S:\user\alexandra\thesis\figures\M2.emf')


%% Fig.M5
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_M5')

figure
set(gcf,'position',[1833         681         476         249])
subplot('position',[0.1 0.1 0.87 0.87])
coord=[];
for i=1:3
    coord=[coord [5 10 26 31]+(i-1)*40];
    if i>1
        line([(i-1)*40,(i-1)*40],[60 140],'linestyle','--','color',[1 1 1]*0.5)
    end
end
hold on
errorbar(coord,peaks,peaks_std,'k.')

axis([0 120 60 125])
set(gca,'xtick',20:40:120,'xticklabel',{'ND4','ND3','ND2'},'ytick',60:20:130,opt_axis,'fontweight','bold')
m=ylabel('Latency, ms');
set(m,opt_label);
line([85 115],[70,70],'color','k')
text(95,73,'30min',opt_text)

saveas(gcf,'S:\user\alexandra\thesis\figures\M5.emf')


%% Fig.M6
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_M6')

figure
set(gcf,'position',[2171         590         597         216])
subplot('position',[0.1 0.1 0.87 0.87])
coord=[];
for i=1:5
    coord=[coord [1:2:24]+(i-1)*29];
    if i>1
        line([(i-1)*29,(i-1)*29],[60 140],'linestyle','--','color',[1 1 1]*0.5)
    end
end
hold on
errorbar(coord,peaks,peaks_std,'k.')

axis([0 145 60 140])
set(gca,'xtick',14.5:29:140,'xticklabel',{'ND6','ND5','ND4','ND2','ND1'},opt_axis,'fontweight','bold')
m=ylabel('Latency, ms');
set(m,opt_label);

line([120 140],[120,120],'color','k')
text(125,123,'20min',opt_text)

saveas(gcf,'S:\user\alexandra\thesis\figures\M6.emf')


%% Fig.M7
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_M7')

figure
set(gcf,'position',[2171         590         597         216])
subplot('position',[0.08 0.17 0.91 0.82])
coord=[];
for i=1:3
    coord=[coord [1:2:24]+(i-1)*29];
    if i>1
        line([(i-1)*29,(i-1)*29],[60 140],'linestyle','--','color',[1 1 1]*0.5)
    end
end
line([i*29,i*29],[60 140],'linestyle','--','color',[1 1 1]*0.5)
coord=[coord 1+3*29];
line([coord(end)+1,coord(end)+1],[60 140],'linestyle','--','color',[1 1 1]*0.5)
for i=1:4
    coord=[coord [1:2:24]+(i-1)*29+coord(37)+1];
    if i>1
        line([(i-1)*29+coord(37)+1,(i-1)*29+coord(37)+1],[60 140],'linestyle','--','color',[1 1 1]*0.5)
    end
end
hold on
errorbar(coord,peaks,peaks_std,'k.')

axis([-3 205 70 140])
set(gca,'xtick',[14.5:29:88 88 103.5:29:200],'xticklabel',{'6','5','4','3','4','3','2','1'},'ytick',70:20:140,opt_axis,'fontweight','bold')
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Latency, ms');
set(m,opt_label);
line([180 200],[120,120],'color','k')
text(182,123,'20min',opt_text)

saveas(gcf,'S:\user\alexandra\thesis\figures\M7.emf')

%% Fig.M8
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_M8')

figure
set(gcf,'position',[2171         590         597         216])
subplot('position',[0.1 0.17 0.89 0.82])
coord=[];
for i=1:5
    coord=[coord [1:9]+(i-1)*13];
    if i>1
        line([(i-1)*13,(i-1)*13],[0 2],'linestyle','--','color',[1 1 1]*0.5)
    end
end
line([i*13,i*13],[0 2],'linestyle','--','color',[1 1 1]*0.5)
coord=[coord 1+5*13];
line([coord(end)+1,coord(end)+1],[0 2],'linestyle','--','color',[1 1 1]*0.5)
for i=1:4
    coord=[coord [1:2:24]+(i-1)*29+coord(46)+1];
    if i>1
        line([(i-1)*29+coord(46)+1,(i-1)*29+coord(46)+1],[0 2],'linestyle','--','color',[1 1 1]*0.5)
    end
end

hold on
errorbar(coord,peaks,peaks_std,'k.')

axis([-2 180 0 1.2])
set(gca,'xtick',[6.5:13:65 66 81.5:29:180],'xticklabel',{'8','7','6','5','4','3','4','3','2','1'},opt_axis,'fontweight','bold')
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Amplitude, normalized');
set(m,opt_label);
line([160 175],[0.5,0.5],'color','k')
text(162,0.55,'15min',opt_text)

saveas(gcf,'S:\user\alexandra\thesis\figures\M8.emf')


%% Fig.M9
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_M9')

figure
set(gcf,'position',[2171         590         597         216])
subplot('position',[0.08 0.17 0.91 0.82])
coord=[];
for i=1:10
    coord=[coord [1:2:24]+(i-1)*29];
    if i>1
        line([(i-1)*29,(i-1)*29],[60 180],'linestyle','--','color',[1 1 1]*0.5)
    end
end
hold on
errorbar(coord,peaks,peaks_std,'k.')

axis([-5 290 60 180])
set(gca,'xtick',14.5:29:280,'xticklabel',{'7','6','5','4','4','4','4','3','2','1'},'ytick',70:20:180,opt_axis,'fontweight','bold')
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Latency, ms');
set(m,opt_label);

line([265 285],[150,150],'color','k')
text(265,155,'20min',opt_text)

saveas(gcf,'S:\user\alexandra\thesis\figures\M9.emf')


%% Fig.K0
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_K0')

figure
set(gcf,'position',[ 2210         386         650         436])
titles=cell(6,1);
titles{1}=sprintf('spontaneous');
titles{2}=sprintf('high contrast');
titles{3}=sprintf('low contrast');
titles{4}=sprintf('spontaneous');
titles{5}=sprintf('high contrast');
titles{6}=sprintf('low contrast');
ax(1:3)=35;
ax(4:6)=25;

corn_posy=[0.56 0.56 0.56 0.08 0.08 0.08];
corn_posx=[0.07 0.39 0.71 0.07 0.39 0.71];
size_x=0.27;
size_y=0.38;
for i=1:6
    subplot(2,3,i,'position',[corn_posx(i),corn_posy(i),size_x,size_y])
    errorbar(data(1,:,i),data(2,:,i),'r')
    hold on
    errorbar(data(3,:,i),data(4,:,i),'b')
    text(4.5,ax(i)*1.05,titles{i},opt_title,'fontsize',14,'Horizontalalignment','center')

    axis([0.5 8.5 0 ax(i)])
    set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},'ytick',0:10:35,opt_axis)
    if i>3
        m=xlabel('ND');
        set(m,opt_label);
    end
    if i==1
        m=ylabel('firing rate mean, Hz');
        set(m,opt_label,'fontsize',12);
    elseif i==4
        m=ylabel('firing rate variance, Hz');
        set(m,opt_label,'fontsize',12);
    end
    m=legend('ON','OFF');
    legend boxoff
    set(m,opt_leg)
end

saveas(gcf,'S:\user\alexandra\thesis\figures\K0.emf')

%% Fig.K01
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_K0')

figure
set(gcf,'position',[ 2137         613         650         223])
titles=cell(2,1);
titles{1}='ON cells';
titles{2}='OFF cells';

a(1,:)=[0.07 0.16 0.44 0.75];
a(2,:)=[0.55 0.16 0.44 0.75];

for i=1:2
    subplot('position',a(i,:))
    errorbar(data((i-1)*2+1,:,1),data((i-1)*2+2,:,1),'k')
    hold on
    errorbar(data((i-1)*2+1,:,2),data((i-1)*2+2,:,2),'r')
    errorbar(data((i-1)*2+1,:,3),data((i-1)*2+2,:,3),'b')
    
    text(4.5,37,titles{i},opt_title,'fontsize',14,'Horizontalalignment','center')
    
    axis([0.5 8.5 0 35])
    set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'}, opt_axis)
    m=xlabel('ND');
    if i==1
        set(m,opt_label);
        m=ylabel('Mean firing rate, Hz');
    end
    set(m,opt_label);
    m=legend({'BG','HC','LC'});
    set(m,opt_leg,'fontsize',10)
    legend boxoff
end

saveas(gcf,'S:\user\alexandra\thesis\figures\K01.emf')


%% Fig.K

clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\summary_all','frMean_spont','onOff','frSTD_spont','frMean_HC',...
    'frSTD_HC','frMean_LC','frSTD_LC','bothContr')

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])
    hold on
    
    b=frMean_HC(onOff>0&bothContr,i)-frMean_LC(onOff>0&bothContr,i);
    [v,m]=sort(frMean_spont(onOff>0&bothContr,i));
    plot(v,b(m),'xr','markersize',5)
    
    b=frMean_HC(onOff<0&bothContr,i)-frMean_LC(onOff<0&bothContr,i);
    [v,m]=sort(frMean_spont(onOff<0&bothContr,i));
    plot(v,b(m),'x','markersize',5)
    
    line([0 250],[0 0],'color','k')
    axis([0 60 -10 30])
    
    text(30,33,['ND',int2str(9-i)],opt_title,'Horizontalalignment','center');

    set(gca,'xtick',[500,1500],'ytick',[500,1500],opt_axis)
    if i<4
        m=ylabel('HC - LC, Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('spontaneous, Hz');
        set(m,opt_label)
    end
    
    set(gca,'xtick',5:10:60,'ytick',-5:10:25, opt_axis)
end

subplot('position',sizes_1)
h=plot(1,1,'xr');
hold on
h1=plot(1,1,'xb');
m=legend({'ON cells','OFF cells'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\K.emf')

%% Fig.K1
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\summary_all','frMean_spont','onOff','frSTD_spont','frMean_HC',...
    'frSTD_HC','frMean_LC','frSTD_LC','bothContr')
figure
set(gcf,'position',[2249         647         433         247])
subplot('position',[0.12 0.15 0.85 0.84])
for i=1:8
    ons=onOff<0&bothContr>0;
    l(i,1)=sum((frMean_HC(ons,i)-frMean_LC(ons,i))<0)/sum(ons);
    ons=onOff>0&bothContr>0;
    l(i,2)=sum((frMean_HC(ons,i)-frMean_LC(ons,i))<0)/sum(ons);
end

hb=bar(l*100);
m=legend(hb,{'OFF','ON'});
legend boxoff
set(m,opt_leg);
axis([0 9 0 70])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},'ytick',10:10:60,opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('% of ON cells');
set(m,opt_label);
saveas(gcf,'S:\user\alexandra\thesis\figures\K1.emf')


%% Fig.K2 - no working version
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\summary_all','frMean_spont','onOff','frSTD_spont','frMean_HC',...
    'frSTD_HC','frMean_LC','frSTD_LC','bothContr')
figure
set(gcf,'position',[1784         566         577         389])
for i=1:8
    subplot(2,4,i)
    offs=onOff<0&bothContr>0;
    [k1,kind1]=hist([frMean_HC(offs,i)-frMean_LC(offs,i); -20; 40],15);
    plot(kind1,k1,'-*')
    off_res(:,i)=frMean_HC(offs,i)-frMean_LC(offs,i);
    hold on
    ons=onOff>0&bothContr>0;
    [k,kind]=hist([frMean_HC(ons,i)-frMean_LC(ons,i); -20; 40],15);
    plot(kind,k,'-*r')
    on_res(:,i)=frMean_HC(ons,i)-frMean_LC(ons,i);
end
for i=1:8
    for j=1:8
        p(i,j)=ranksum(on_res(:,i),on_res(:,j));
    end
end
figure
set(gcf,'position',[1784         566         577         389])
for i=1:8
    subplot(2,4,i)
    offs=onOff<0&bothContr>0;
    [k1,kind1]=hist([frMean_HC(offs,i)-frMean_LC(offs,i); -20; 40],15);
    plot(kind1,k1/sum(offs)*100,'-*')
    off_res(:,i)=frMean_HC(offs,i)-frMean_LC(offs,i);
    hold on
    ons=onOff>0&bothContr>0;
    [k,kind]=hist([frMean_HC(ons,i)-frMean_LC(ons,i); -20; 40],15);
    plot(kind,k/sum(ons)*100,'-*r')
    on_res(:,i)=frMean_HC(ons,i)-frMean_LC(ons,i);
end

figure
off_res=cell(8,1);
on_res=cell(8,1);
for i=1:8
    subplot(2,4,i)
    offs=onOff<0&bothContr>0;
    a=abs(frMean_LC(offs,i)./frMean_HC(offs,i));
    a(isnan(a))=[];
    a(isinf(a))=[];
    a(a>20)=[];
   off_res{i}=a;
    plot(a,'o')
    hold on
    ons=onOff>0&bothContr>0;
    a=abs(frMean_LC(ons,i)./frMean_HC(ons,i));
        a(isnan(a))=[];
    a(isinf(a))=[];
    a(a>20)=[];
    on_res{i}=a;
    plot(a,'or')
end
for i=1:8
    for j=1:8
        p(i,j)=ranksum(on_res{i},on_res{j});
    end
end

%% Fig.E
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_E')

figure
set(gcf,'position',[1784         201        1101         754])

for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,1,i),'r');
    hold on
    plot(data(:,2,i),'b');
    
    plot(4701:9200,data(:,3,i),'r');
    plot(4701:9200,data(:,4,i),'b');

    axis([0 9200 -25 65])

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2.5','0.5','2.5'},opt_axis)
    
    if i==1||i==2
        m=ylabel('firing rate, Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end

    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55],'color','k')
    line([500 2500],[60,60],'color','k')
    line([2500 4500],[55,55],'color','k')
    line([500 500],[55,60],'color','k')
    line([2500 2500],[55,60],'color','k')
    
    line([1 500]+4700,[55,55],'color','k')
    line([500 2500]+4700,[50,50],'color','k')
    line([2500 4500]+4700,[55,55],'color','k')
    line([500 500]+4700,[50,55],'color','k')
    line([2500 2500]+4700,[50,55],'color','k')
end
subplot(3,3,1,'position',sizes_1)
h=plot(1,1,'r');
hold on
h1=plot(1,1,'b');
m=legend({'decreasing','increasing'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\E.emf')


%% Fig.E2 
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\chirp_for_E2.mat','chirp_acc')

figure
set(gcf,'position',[2210         281         650         541])

corn_posy=[0.55 0.31 0.07 0.79 0.55 0.31 0.07];
corn_posx=[0.06 0.06 0.06 0.54 0.54 0.54 0.54]+0.01;
size_x=0.43;
size_y=0.17;
sizes_1=[0.06 0.76 0.44 0.23];

ax=cell(1,14);
cc=1;
for i=2:2:25
    ax{1,cc}=int2str(i);
    cc=cc+1;
end

for i=2:8
    subplot('position',[corn_posx(i-1),corn_posy(i-1),size_x,size_y])

    plot(chirp_acc(:,1,i),'b');
    hold on    
    plot(chirp_acc(:,2,i),'r');

    line([0 22000],[0 0],'color','k')
    axis([2500 21000 -30 40])

    text(12000,45,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    set(gca,'xtick',2000:2000:25000,'xticklabel',ax,opt_axis)
    if i==4||i==8
        m=xlabel('time,s');
        set(m,opt_label);
    end
    if i<5
        m=ylabel('firing rate');
        set(m,opt_label);
    end    
end

subplot('position',sizes_1)
h=plot(1,1,'b');
hold on
h1=plot(1,1,'r');
m=legend({sprintf('typical\nON cells'),sprintf('non-typical\nON cells')});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\E2.emf')


%% Fig.E3 
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\chirp_for_E3.mat','chirp_acc')

figure
set(gcf,'position',[2210         281         650         541])

corn_posy=[0.55 0.31 0.07 0.79 0.55 0.31 0.07];
corn_posx=[0.06 0.06 0.06 0.54 0.54 0.54 0.54]+0.01;
size_x=0.43;
size_y=0.17;
sizes_1=[0.06 0.76 0.44 0.23];

ax=cell(1,14);
cc=1;
for i=2:2:25
    ax{1,cc}=int2str(i);
    cc=cc+1;
end

for i=2:8
    subplot('position',[corn_posx(i-1),corn_posy(i-1),size_x,size_y])

    plot(chirp_acc(:,i),'b');

    line([0 22000],[0 0],'color','k')
    axis([2500 21000 -10 65])

    text(12000,72,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    set(gca,'xtick',2000:2000:25000,'xticklabel',ax,opt_axis)
    if i==4||i==8
        m=xlabel('time,s');
        set(m,opt_label);
    end
    if i<5
        m=ylabel('firing rate');
        set(m,opt_label);
    end    
end

saveas(gcf,'S:\user\alexandra\thesis\figures\E3.emf')


%% Fig.N
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_N.mat')

figure
set(gcf,'position',[2186         156         533         686])
nds='666555444333222';

corn_posy=[ 0.81 0.81 0.81 0.62 0.62 0.62 0.43 0.43 0.43 0.24 0.24 0.24 0.05 0.05 0.05]+0.01;
corn_posx=[0.04 0.36 0.68 0.04 0.36 0.68  0.04 0.36 0.68 0.04 0.36 0.68  0.04 0.36 0.68 ]+0.01;
size_x=0.29;
size_y=0.14;
sizes_1=[0.06 0.72 0.287 0.29];
cm=3;
for i=1:3:15
    subplot(5,3,i,'position',[corn_posx(i),corn_posy(i),size_x,size_y])

    plot(data_acc(:,cm,1),'b','linewidth',2)
    hold on
    plot(data_acc(:,cm+1,1),'r','linewidth',2)

    text(150,1.7,['ND',nds(i)],opt_title,'horizontalalignment','center');

    line([0 300],[0 0],'color','k')
    axis([0 300 -1.5 1.5])
    set(gca,'xtick',50:100:250,'xticklabel',{'50','150','250'},opt_axis)
    set(gca,'ytick',0,'yticklabel','',opt_axis)
    if i>12
        m=xlabel('time, ms');
        set(m,opt_label);
    end
    m=ylabel('a.u');
    set(m,opt_label);
    
    subplot(5,3,i+1,'position',[corn_posx(i+1),corn_posy(i+1),size_x,size_y])
    plot(data_acc(:,cm,2),'b','linewidth',2)
    hold on
    plot(data_acc(:,cm+1,2),'r','linewidth',2)

    text(150,1.7,['ND',nds(i)],opt_title,'horizontalalignment','center');
    
    line([0 300],[0 0],'color','k')
    axis([0 300 -1.5 1.5])
    if i==4
        axis([0 300 -2.5 1.5])
    end
    set(gca,'xtick',50:100:250,'xticklabel',{'50','150','250'},opt_axis)
    set(gca,'ytick',0,'yticklabel','',opt_axis)
    if i>12
        m=xlabel('time, ms');
        set(m,opt_label);
    end

    subplot(5,3,i+2,'position',[corn_posx(i+2),corn_posy(i+2),size_x,size_y])
    plot(data_acc(:,cm,3),'b','linewidth',2)
    hold on
    plot(data_acc(:,cm+1,3),'r','linewidth',2)

    text(150,1.7,['ND',nds(i)],opt_title,'horizontalalignment','center');
    
    line([0 300],[0 0],'color','k')
    axis([0 300 -1.5 1.5])
   
    set(gca,'xtick',50:100:250,'xticklabel',{'50','150','250'},opt_axis)
    set(gca,'ytick',0,'yticklabel','',opt_axis)
    if i>12
        m=xlabel('time, ms');
        set(m,opt_label);
    end
    i
    cm=cm+2;
end

saveas(gcf,'S:\user\alexandra\thesis\figures\N.emf')

%% Fig.N_alt_alt
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_N_alt_alt.mat')

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data_acc_on(:,2,i),data_acc_on(:,1,i),'-*r')
    hold on
    plot(data_acc_off(:,2,i),data_acc_off(:,1,i),'-*b')

    axis([-40 20 0 40]) 
    line([0,0],[0 40],'color','k')
    
    set(gca,'xtick',-40:10:20,opt_axis)
    
    text(-10,42,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    if i<4
        m=ylabel('% of cells');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('Lag');
        set(m,opt_label)
    end

end

subplot('position',sizes_1)
h=plot(1,1,'r');
hold on
h1=plot(1,1,'b');
m=legend({'ON cells','OFF cells'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\N_alt_alt.emf')

%% Fig.N4
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_N4.mat')

figure
set(gcf,'position',[2249         647         433         247])
subplot('position',[0.12 0.15 0.85 0.84])

errorbar(ons(:,1),ons(:,2),'r')
hold on
errorbar(offs(:,1),offs(:,2))

m=legend('ON cells','OFF cells','location','best');
set(m,opt_leg);
legend('boxoff')

tt=max([ons(:,1)';offs(:,1)']);
tmp=1:8;tmp(p>=0.05)=[];
plot(tmp,tt(tmp)+1,'*k')

axis([0 9 2 11])
line([0 9],[5,5],'color','k')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('gain ratio');
set(m,opt_label);


saveas(gcf,'S:\user\alexandra\thesis\figures\N4.emf')


%% Fig.N6
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_N6.mat')

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data_acc(:,1,i),data_acc(:,2,i),'-*r')
    hold on
    plot(data_acc(:,3,i),data_acc(:,4,i),'-*')
    line([5,5],[0 50],'color','k')
    
    axis([0 20 0 50])
    
    set(gca,'xtick',0:5:20,opt_axis)
    
    text(10,54,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    if i<4
        m=ylabel('% of cells');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('Gain ratio');
        set(m,opt_label)
    end

end

subplot('position',sizes_1)
h=plot(1,1,'r');
hold on
h1=plot(1,1,'b');
m=legend({'ON cells','OFF cells'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\N6.emf')


%% Fig.N5
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_N5.mat')

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data{1,i}(:,1),data{1,i}(:,2),'xr')
    hold on
    plot(data{2,i}(:,1),data{2,i}(:,2),'xb')

    line([-30 30],[0 0],'color','k')
    line([0,0],[-30 30],'color','k')
    axis([-30 30 -30 30])
    
    set(gca,'xtick',-40:10:20,opt_axis)
    
    text(0,35,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    if i<4
        m=ylabel('FR change,HZ');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('Lag,ms');
        set(m,opt_label)
    end

end

subplot('position',sizes_1)
h=plot(1,1,'r');
hold on
h1=plot(1,1,'b');
m=legend({'ON cells','OFF cells'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\N5.emf')

%% Fig.N7
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_N7.mat')


corn_posy=[0.35    0.0200    0.68    0.35    0.0200    0.68    0.35    0.0200];
size_y=0.27;

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    errorbar(1,res(1,i,1),res(1,i,2),'.k')
    hold on
    errorbar(2,res(2,i,1),res(2,i,2),'.k')
    bar(res(1,i,1),'facecolor','b')
    bar(2,res(2,i,1),'facecolor','r')
    
    axis([0 3 -8 11])
    set(gca,'xtick',0,'xticklabel','',opt_axis)
    
    text(1.5,12.3,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    if i<4
        m=ylabel('Relative input');
        set(m,opt_label)
    end

end

subplot('position',sizes_1)
h=bar(1,1,'facecolor','b')
hold on
h1=bar(1,1,'facecolor','r')
m=legend({sprintf('typical\nON cells'),sprintf('non-typical\nON cells')});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\N7.emf')

%% Fig.N8
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\chirp_for_E2.mat','chirp_acc')
a=chirp_acc; % typical and non-typical ON cells
load('S:\data\alexandra\MEA_data\analysis\chirp_for_E3.mat','chirp_acc') % off cells

for i=1:8
    typOff(i)=corr(a(:,1,i),chirp_acc(:,i)); % typical and off cells correlation
    nontypOff(i)=corr(a(:,2,i),chirp_acc(:,i)); % non-typical and off cells correlation
end

figure
set(gcf,'position',[2249         609         433         285])
subplot('position',[0.12 0.15 0.85 0.83])

bar([typOff; nontypOff]')
axis([0 9 -1 0])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis,'fontweight','bold')

m=xlabel('ND');
set(m,opt_label);
m=ylabel('Correlation');
set(m,opt_label);
m=legend({'typical ON vs OFF', 'non-typical ON vs OFF'},'location','southeast');
set(m,opt_leg)
legend boxoff

saveas(gcf,'S:\user\alexandra\thesis\figures\N8.emf')

%% Fig.N9
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_N9.mat')

bord=2; % threshold to cut typical and non-typical cells
for i=1:8
    a=data{1,i}(:,2); % firing rate change, on cells
    b=data{1,i}(:,1); % delay, on cells
    
    anontyp=a<-bord;
    res(1,1,i)=mean(b(anontyp));
    res(1,2,i)=std(b(anontyp))/sqrt(sum(anontyp));
    [~, p(1,i)]=ttest(b(anontyp)); 
    
    atyp=a>bord;
    res(2,1,i)=mean(b(atyp));
    res(2,2,i)=std(b(atyp))/sqrt(sum(atyp));
    [~, p(2,i)]=ttest(b(atyp));   
    t=b(atyp);
        
    a=data{2,i}(:,2); % firing rate change, off cells
    b=data{2,i}(:,1); % delay, off cells
    atyp=a>bord;
    res(3,1,i)=mean(b(atyp));
    res(3,2,i)=std(b(atyp))/sqrt(sum(atyp));
    [~, p(3,i)]=ttest(b(atyp));   
    
    j(i)=ranksum(t,b(atyp));
end

corn_posy=[0.35    0.0200    0.68    0.35    0.0200    0.68    0.35    0.0200];
size_y=0.27;
figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    errorbar(res(:,1,i),res(:,2,i),'.k')
    hold on
    bar(1,res(1,1,i),'facecolor',[0.3 0.8 0.1]); % first plot ON non-typical
    bar(2,res(2,1,i),'facecolor','r'); % plot ON typical
    bar(3,res(3,1,i),'facecolor','b'); % plot OFF
    
    for j=1:3
        if p(j,i)<0.01
            if res(j,1,i)<0
                plot(j-0.05,res(j,1,i)-res(j,2,i)-1,'*k')
                plot(j+0.05,res(j,1,i)-res(j,2,i)-1,'*k')
            else
                plot(j-0.05,res(j,1,i)+res(j,2,i)+1,'*k')
                plot(j+0.05,res(j,1,i)+res(j,2,i)+1,'*k')
            end
        end
    end
    
    axis([0 4 -14 6])
    
    set(gca,'xtick',0,'xticklabel','',opt_axis)
    
    text(2,7.5,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    if i<4
        m=ylabel('Lag,ms');
        set(m,opt_label)
    end
end

subplot(3,3,1,'position',sizes_1)
h=bar(1,1,'facecolor',[0.3 0.8 0.1]);
hold on
h1=bar(1,1,'facecolor','r');
h2=bar(1,1,'facecolor','b');
m=legend({'non-typical ON','typical ON','OFF'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
delete(h2)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\N9.emf')


%% Fig.X
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_X.mat')

figure
set(gcf,'position',[1784         201        1101         754])
k=70;
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(4701:9200,data(:,1,i),'k');
    hold on
    hold on
    plot(data(:,2,i),'k');

    axis([0 9200 -5 140])

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2.5','0.5','2.5'},opt_axis)
    if i==1||i==2
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end

    
    line([1 500],[55,55]+k,'color','k')
    line([500 2500],[60,60]+k,'color','k')
    line([2500 4500],[55,55]+k,'color','k')
    line([500 500],[55,60]+k,'color','k')
    line([2500 2500],[55,60]+k,'color','k')
    
    line([1 500]+4700,[55,55]+k,'color','k')
    line([500 2500]+4700,[50,50]+k,'color','k')
    line([2500 4500]+4700,[55,55]+k,'color','k')
    line([500 500]+4700,[50,55]+k,'color','k')
    line([2500 2500]+4700,[50,55]+k,'color','k')
end


saveas(gcf,'S:\user\alexandra\thesis\figures\X.emf')

%% Fig.Y
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Y.mat')

figure
set(gcf,'position',[1784         201        1101         754])
k=85;
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(4701:9200,data(:,1,i),'k');
    hold on
    plot(data(:,2,i),'k');

    axis([0 9200 -5 150])

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2.5','0.5','2.5'},opt_axis)
    if i==1||i==2
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
    
    line([1 500],[55,55]+k,'color','k')
    line([500 2500],[60,60]+k,'color','k')
    line([2500 4500],[55,55]+k,'color','k')
    line([500 500],[55,60]+k,'color','k')
    line([2500 2500],[55,60]+k,'color','k')
    
    line([1 500]+4700,[55,55]+k,'color','k')
    line([500 2500]+4700,[50,50]+k,'color','k')
    line([2500 4500]+4700,[55,55]+k,'color','k')
    line([500 500]+4700,[50,55]+k,'color','k')
    line([2500 2500]+4700,[50,55]+k,'color','k')
end


saveas(gcf,'S:\user\alexandra\thesis\figures\Y.emf')

%% Fig.Z
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z.mat')

figure
set(gcf,'position',[1784         201        1101         754])
k=10;
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,1,i),'k');
    hold on
    hold on
    plot(4701:9200,data(:,2,i),'k');

    axis([0 9200 -5 80])

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2.5','0.5','2.5'},opt_axis)
    if i==1||i==2
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
    
    line([1 500],[55,55]+k,'color','k')
    line([500 2500],[60,60]+k,'color','k')
    line([2500 4500],[55,55]+k,'color','k')
    line([500 500],[55,60]+k,'color','k')
    line([2500 2500],[55,60]+k,'color','k')
    
    line([1 500]+4700,[55,55]+k,'color','k')
    line([500 2500]+4700,[50,50]+k,'color','k')
    line([2500 4500]+4700,[55,55]+k,'color','k')
    line([500 500]+4700,[50,55]+k,'color','k')
    line([2500 2500]+4700,[50,55]+k,'color','k')
end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z.emf')

%% Fig.Z2_alt
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z2_alt.mat')

figure
set(gcf,'position',[2139         606         636         179])
corn_posy=0.2;
corn_posx=[0.08 0.39 0.7];
size_x=0.26;
size_y=0.65;
for i=1:3
    subplot(1,3,i,'position',[corn_posx(i),corn_posy,size_x,size_y])

    plot(data(:,i),'linewidth',2,'color','k')
    axis([0 4500 0 45])   
    line([1 500],[90,90]-50,'color','k')
    line([500 2500],[85,85]-50,'color','k')
    line([2500 4500],[90,90]-50,'color','k')
    line([500 500],[85,90]-50,'color','k')
    line([2500 2500],[85,90]-50,'color','k')
    set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'},opt_axis)
    if i==1
        m=ylabel('Hz');
        set(m,opt_label)
    end
    m=xlabel('time,s');
    set(m,opt_label)
    m=title(['ND',int2str(7-i)]);
    set(m,opt_title);
end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z2_alt.emf')

%% Fig.Z1_alt
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z1_alt.mat')

figure
set(gcf,'position',[2139         606         636         179])
corn_posy=0.2;
corn_posx=[0.08 0.39 0.7];
size_x=0.26;
size_y=0.65;
ord=[2 3 1];
titles=cell(3,1);
titles{1}='no ON response';
titles{2}='early ON response';
titles{3}='delayed ON response';
for i=1:3
    subplot(1,3,i,'position',[corn_posx(i),corn_posy,size_x,size_y])

    plot(data(:,ord(i)),'linewidth',2,'color','k')
    line([1 500],[140,140],'color','k')
    line([500 2500],[130,130],'color','k')
    line([2500 4500],[140,140],'color','k')
    line([500 500],[130,140],'color','k')
    line([2500 2500],[140,130],'color','k')
    axis([0 4500 0 150])
    set(gca,'xtick',500:1000:3500,'xticklabel',{'0.5','1.5','2.5','3.5'},opt_axis)
    if i==1
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    m=xlabel('time,s');
    set(m,opt_label)
    m=title(titles{i});
    set(m,opt_title);
end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z1_alt.emf')


%% Fig.Z4
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z4.mat')
figure
set(gcf,'position',[2357         706         484         217])
subplot('position',[0.1 0.15 0.89 0.85])
bar(data(1:3))
hold on
bar(6:8,data(4:6))
axis([0 9 0 70])
set(gca,'xtick',[1:3 6:8],'xticklabel',{'ND6','ND5','ND4','ND6','ND5','ND4'},'ytick',10:20:60,opt_axis)
m=text(0.8,-10,'Early ON response',opt_text);
m=text(5.6,-10,'Delayed ON response',opt_text);
m=ylabel('% of OFF cells');
set(m,opt_label);

saveas(gcf,'S:\user\alexandra\thesis\figures\Z4.emf')

%% Fig.Z5
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z5.mat')
figure
set(gcf,'position',[2357         706         401         217])
subplot('position',[0.15 0.2 0.84 0.79])
bar([1 2.5],data,'BarLayout','stacked')
text(1,-10,sprintf('Early\nON response'),'HorizontalAlignment','center',opt_text)
text(2.5,-10,sprintf('Delayed\nON response'),'HorizontalAlignment','center',opt_text)
set(gca,'xtick',0,'xticklabel','','ytick',20:20:60,opt_axis,'fontweight','bold')
m=ylabel('% of OFF cells');
set(m,opt_label);
axis([0 5.5 0 70])
m=legend(sprintf('always\npresent'),sprintf('never\npresent'));
set(m,opt_leg);
legend boxoff

saveas(gcf,'S:\user\alexandra\thesis\figures\Z5.emf')

%% Fig.Z8_alt
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z8_alt.mat')


corn_posx=[0.0800    0.0800    0.39    0.39    0.39    0.70   0.70    0.70];
size_x=0.27;
figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,1,i),data(:,2,i),'-*')
    hold on
    plot(data(:,3,i),data(:,4,i),'-*r')   
    
    axis([-1 1 0 80])
    
    text(0,90,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

     set(gca,'xtick',[-0.5 0 0.5],opt_axis)
    
    if i<4
        m=ylabel('% of OFF cells');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('correlation coefficient');
        set(m,opt_label)
    end

end
sizes_1=[0.0600    0.7200    0.270    0.2500];
subplot(3,3,1,'position',sizes_1)
h=plot(1,1,'r');
hold on
h1=plot(1,1,'b');
m=legend({'ON response','OFF response'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\Z8_alt.emf')

%% Fig.Z6_alt
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z6_alt.mat')

figure
set(gcf,'position',[1784         201        1101         754])
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])
    
    plot(data{i}(1,:),data{i}(2,:),'o')
    axis([-50 50 -0 2])
    
    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[-25,0,25],'xticklabel',{'-25','0','25'},opt_axis)
    if i==1||i==2
        m=ylabel('amplitude ratio');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('latency delta, ms');
        set(m,opt_label)
    end

    line([0,0],[0,2],'color','k')
    line([-50 50],[1,1],'color','k')
end

subplot(3,3,1,'position',sizes_1)
h1=plot(1,1,'ob');
m=legend({'OFF cells'});
set(m,opt_title)
legend('boxoff')
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\Z6_alt.emf')


%% Fig.Z14
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z6_alt.mat')

for i=1:8
    % cut out noise
%     cc=abs(data{i}(1,:))<55&data{i}(2,:)<2&data{i}(2,:)>0.2;
    cc=abs(data{i}(1,:))<25&data{i}(2,:)<1.5&data{i}(2,:)>0.5;
    
    % amplitudes ratio
    tmp=data{i}(2,cc);
    k(i)=nanmean(tmp);
    ks(i)=nanstd(tmp)/sqrt(size(tmp,2));
    [n(i),b(i)]=ttest(tmp-1);
    
    r(i)=kstest(tmp);
    
    % latency difference
    tmp1=data{i}(1,cc);
    t(i)=nanmean(tmp1);
    ts(i)=nanstd(tmp1)/sqrt(size(tmp1,2));
    tmptmp{i}=tmp1;
    [m(i),a(i)]=ttest(tmp1);  

    r1(i)=kstest(tmp1);
    
    % correlation of latency and amplitude
    kp(i)=corr(tmp',tmp1');
end

for i=1:8
    for j=1:8
        ll(i,j)=ranksum(tmptmp{i},tmptmp{j})
    end
end




figure
set(gcf,'position',[ 2137         613         650         223])

subplot('position',[0.07 0.16 0.41 0.83])
errorbar(k,ks,'k')
hold on
line([0 9],[1,1],'color','k')
for i=1:8
    if b(i)<0.001
        plot(i-0.1,k(i)+ks(i)+0.05,'*k')
        plot(i+0.1,k(i)+ks(i)+0.05,'*k')
    end
end
axis([0 9 0.7 1.3])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Amplitude ratio');
set(m,opt_label);

subplot('position',[0.57 0.16 0.41 0.83])
errorbar(t,ts,'k')
hold on
for i=1:8
    if a(i)<0.001
        plot(i-0.1,t(i)+ts(i)+0.5,'*k')
        plot(i+0.1,t(i)+ts(i)+0.5,'*k')
    elseif  a(i)<0.05
        plot(i,t(i)+ts(i)+0.5,'*k')     
    end
end
line([0 9],[0,0],'color','k')
axis([0 9 -13 3])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Latency difference,ms');
set(m,opt_label);

saveas(gcf,'S:\user\alexandra\thesis\figures\Z14.emf')


%% Fig.Z7
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z7.mat')

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,i,1),data(:,i,2),'o','markersize',2)
    
    axis([0 2000 0 2000])
    
    text(850,2200,['ND',int2str(9-i)],opt_title);
%     m=title(['ND',int2str(9-i)]);
%     set(m,opt_title);
    set(gca,'xtick',[500,1500],'ytick',[500,1500],opt_axis)
    if i<4
        m=ylabel('Transiency White');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('Transiency Black');
        set(m,opt_label)
    end
    line([0 2000],[0 2000],'color','k')
    p(i)=ranksum(data(:,i,1),data(:,i,2));
end


saveas(gcf,'S:\user\alexandra\thesis\figures\Z7.emf')

%% Fig.Z9
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z9.mat')

figure
set(gcf,'position',[1784         201        1101         754])
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])


    plot(data(:,i,1),'color','k','linewidth',2)    
    hold on
    plot(4701:9200,data(:,i,2),'color','k','linewidth',2)
    
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

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2','0.5','2'},opt_axis)
    
    if i==1||i==2
        m=ylabel('firing rate, Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z9.emf')


%% Fig.Z4_comb

clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z4.mat')
a=data;
load('S:\data\alexandra\MEA_data\analysis\data_for_Z4_and_Z5_white.mat')
figure
set(gcf,'position',[2357         598         520         325])
subplot('position',[0.1 0.15 0.85 0.8])
bar([a(1:3); data(1:3)]')
hold on
bar(6:8,[a(4:6); data(4:6)]')
axis([0 9 0 70])
m=legend({'black step','white step'},'location','northwest');
set(m,opt_leg);
set(gca,'xtick',[1:3 6:8],'xticklabel',{'ND6','ND5','ND4','ND6','ND5','ND4'},'ytick',20:20:70,opt_axis,'fontweight','bold')
m=text(0.5,-8,'Early ON response');
set(m,'fontsize',14,'fontweight','bold','fontname','arial')
m=text(5.5,-8,'Delayed ON response');
set(m,'fontsize',14,'fontweight','bold','fontname','arial')
m=ylabel('% of OFF cells');
set(m,opt_label);
legend boxoff

saveas(gcf,'S:\user\alexandra\thesis\figures\Z4_comb.emf')


%% Fig.Z5_comb

clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z5.mat')
a=data;
load('S:\data\alexandra\MEA_data\analysis\data_for_Z4_and_Z5_white.mat')
b=data1;

figure
set(gcf,'position',[2357         601         430         322])
subplot('position',[0.1 0.15 0.89 0.84])

bar([1,3.3],a,0.2,'BarLayout','stacked')
hold on
bar([1.8,4.1],b,0.2,'BarLayout','stacked')
set(gca,'xtick',[1 1.8 3.3 4.1],'xticklabel',{'black step','white step','black step','white step'},'ytick',20:20:70,opt_axis,'fontweight','bold')

m=text(0.65,-8,'Early ON response');
set(m,'fontsize',12,'fontweight','bold','fontname','arial')
m=text(2.7,-8,'Delayed ON response');
set(m,'fontsize',12,'fontweight','bold','fontname','arial')

m=ylabel('% of OFF cells');
set(m,opt_label);
axis([0.5 4.7 0 70])
m=legend('always present','never present');
legend boxoff
set(m,opt_leg);

saveas(gcf,'S:\user\alexandra\thesis\figures\Z5_comb.emf')



%% Fig.Z10
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z10.mat')

figure
set(gcf,'position',[1784         201        1101         754])
k=5;
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])


    plot(data(:,i,1),'color','k','linewidth',2)    
    hold on
    plot(4701:9200,data(:,i,2),'color','k','linewidth',2)
    
    axis([1 9200 0 75])
   
    line([1 500],[55,55]+k,'color','k')
    line([500 2500],[60,60]+k,'color','k')
    line([2500 4500],[55,55]+k,'color','k')
    line([500 500],[55,60]+k,'color','k')
    line([2500 2500],[55,60]+k,'color','k')
    
    line([1 500]+4700,[55,55]+k,'color','k')
    line([500 2500]+4700,[50,50]+k,'color','k')
    line([2500 4500]+4700,[55,55]+k,'color','k')
    line([500 500]+4700,[50,55]+k,'color','k')
    line([2500 2500]+4700,[50,55]+k,'color','k')

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2','0.5','2'},opt_axis)
    
    if i==1||i==2
        m=ylabel('firing rate, Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z10.emf')


%% Fig.Z15
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z15.mat')

figure
set(gcf,'position',[1784         201        1101         754])
k=35;
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])


    plot(data(:,i,1),'color','k','linewidth',2)    
    hold on
    plot(4701:9200,data(:,i,2),'color','k','linewidth',2)
    
    axis([1 9200 0 100])
   
    line([1 500],[55,55]+k,'color','k')
    line([500 2500],[60,60]+k,'color','k')
    line([2500 4500],[55,55]+k,'color','k')
    line([500 500],[55,60]+k,'color','k')
    line([2500 2500],[55,60]+k,'color','k')
    
    line([1 500]+4700,[55,55]+k,'color','k')
    line([500 2500]+4700,[50,50]+k,'color','k')
    line([2500 4500]+4700,[55,55]+k,'color','k')
    line([500 500]+4700,[50,55]+k,'color','k')
    line([2500 2500]+4700,[50,55]+k,'color','k')

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2','0.5','2'},opt_axis)
    
    if i==1||i==2
        m=ylabel('firing rate, Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z15.emf')


%% Fig.Z11
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z11.mat')

figure
set(gcf,'position',[1784         201        1101         754])

for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])


    plot(lat(:,i),ampl(:,i),'xk')    
    
    axis([1 500 0 150])
   
    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',0:200:400,opt_axis)
    
    if i==1||i==2
        m=ylabel('amplitude,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('latency, ms');
        set(m,opt_label)
    end
end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z11.emf')


%% Fig.Z12
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z11.mat')

figure
set(gcf,'position',[ 2137         613         650         223])

subplot('position',[0.07 0.16 0.41 0.75])
errorbar(mean(ampl),std(ampl)/sqrt(200),'k')
axis([0 9 0 70])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis,'fontweight','bold')
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Amplitude,Hz');
set(m,opt_label);

subplot('position',[0.57 0.16 0.41 0.75])
errorbar(mean(lat),std(lat)/sqrt(200),'k')
axis([0 9 140 270])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis,'fontweight','bold')
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Latency, ms');
set(m,opt_label);

saveas(gcf,'S:\user\alexandra\thesis\figures\Z12.emf')


%% Fig.Z13
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z13.mat')

corn_posx=[0.0800    0.0800    0.39    0.39    0.39    0.70   0.70    0.70];
size_x=0.27;
figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,1,i),data(:,2,i),'-*k')    
    
    axis([1 2000 0 50])
    
    text(1000,55,['ND',int2str(9-i)],opt_title,'horizontalalignment','center');

    set(gca,'xtick',0:500:1500,opt_axis)
    
    if i<4
        m=ylabel('% of OFF cells');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('transiency,ms');
        set(m,opt_label)
    end

end

saveas(gcf,'S:\user\alexandra\thesis\figures\Z13.emf')


%% Fig.R
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R.mat')

figure
set(gcf,'position',[1784         201        1101         754])
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,1,i),'k');
    hold on
    hold on
    plot(4701:9200,data(:,2,i),'k');

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2.5','0.5','2.5'},opt_axis)
    if i==1||i==2
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
    
    axis([1 9200 -20 50])
    
    line([1 500],[40,40],'color','k')
    line([500 2500],[45,45],'color','k')
    line([2500 4500],[40,40],'color','k')
    line([500 500],[40,45],'color','k')
    line([2500 2500],[40,45],'color','k')
    
    line([1 500]+4700,[40,40],'color','k')
    line([500 2500]+4700,[35,35],'color','k')
    line([2500 4500]+4700,[40,40],'color','k')
    line([500 500]+4700,[35,40],'color','k')
    line([2500 2500]+4700,[35,40],'color','k')
end

saveas(gcf,'S:\user\alexandra\thesis\figures\R.emf')

%% Fig.R1
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R1.mat')

figure
set(gcf,'position',[ 1870         542         935         243])
corn_posy=0.2;
corn_posx=[0.08 0.39 0.7];
size_x=0.26;
size_y=0.65;
for i=1:3
    subplot(1,3,i,'position',[corn_posx(i),corn_posy,size_x,size_y])

    plot(data(:,i),'linewidth',2,'color','k')
    axis([0 4500 0 40])   
    line([1 500],[32,32],'color','k')
    line([500 2500],[37,37],'color','k')
    line([2500 4500],[32,32],'color','k')
    line([500 500],[32,37],'color','k')
    line([2500 2500],[32,37],'color','k')
    set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'},opt_axis)
    if i==1
        m=ylabel('Hz');
        set(m,opt_label)
    end
    m=xlabel('time,s');
    set(m,opt_label)
    m=title(['ND',int2str(7-i)]);
    set(m,opt_title);
end

saveas(gcf,'S:\user\alexandra\thesis\figures\R1.emf')

%% Fig.R2
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R2.mat')

figure
set(gcf,'position',[ 1870         542         935         243])
corn_posy=0.2;
corn_posx=[0.08 0.39 0.7];
size_x=0.26;
size_y=0.7;
for i=1:3
    subplot(1,3,i,'position',[corn_posx(i),corn_posy,size_x,size_y])

    plot(data(:,i),'linewidth',2,'color','k')
    line([1 500],[140,140],'color','k')
    line([500 2500],[130,130]+20,'color','k')
    line([2500 4500],[140,140],'color','k')
    line([500 500],[130,140]+10,'color','k')
    line([2500 2500],[140,130]+10,'color','k')
    axis([0 4500 0 160])
    set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'},opt_axis)
    if i==1
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    m=xlabel('time,s');
    set(m,opt_label)
end

saveas(gcf,'S:\user\alexandra\thesis\figures\R2.emf')

%% Fig.R3
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R3_R4_R5_R6.mat','data','bardata')

figure
set(gcf,'position',[1784         201        1101         754])
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])
    
    hold on
    plot(data(:,1,i),'color',[0.3 0.8 0.1],'linewidth',2)
    plot(data(:,2,i),'linewidth',2)
    plot(data(:,3,i),'r','linewidth',2)
    bar(4650,bardata(1,i),180,'facecolor',[0.3 0.8 0.1])
    bar(4850,bardata(2,i),180,'facecolor','b')
    bar(5050,bardata(3,i),180,'facecolor','r')

    axis([0 5250 -15 55])
    
    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500],'xticklabel',{'0.5','2.5'},opt_axis)
    if i==1||i==2
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
    line([1 500],[50,50]-5,'color','k')
    line([500 2500],[55,55]-5,'color','k')
    line([2500 4500],[50,50]-5,'color','k')
    line([500 500],[55,50]-5,'color','k')
    line([2500 2500],[55,50]-5,'color','k')
end
subplot(3,3,1,'position',sizes_1)
h=plot(1,1,'color',[0.3 0.8 0.1]);
hold on
h1=plot(1,1,'b');
h2=plot(1,1,'r');
m=legend({'transient','gap-like','sustained'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
delete(h2)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\R3.emf')


%% Fig.R3a
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R3_R4_R5_R6.mat','data','bardata')

figure
set(gcf,'position',[2162         522         606         293])
subplot('position',[0.08 0.10 0.91 0.89])
a=bar(bardata','hist');
set(a(1),'facecolor',[0.3 0.8 0.1])
set(a(2),'facecolor','b')
set(a(3),'facecolor','r')
set(gca,'xtick',1:8,'xticklabel',{'ND8','ND7','ND6','ND5','ND4','ND3','ND2','ND1'},'ytick',10:20:70,opt_axis)
m=ylabel('%');
set(m,opt_label);
axis([0.3 8.7 0 75])
m=legend({'trans','gap','sust'},'location','northwest');
set(m,opt_leg)
legend boxoff

saveas(gcf,'S:\user\alexandra\thesis\figures\R3a.emf')

%% Fig.R4
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R3_R4_R5_R6.mat','spont','spont_std','p')

figure
set(gcf,'position',[ 2210         386         650         436])

for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    errorbar(100:100:300,spont(i,:),spont_std(i,:),'.k')
    hold on
    bar(100,spont(i,1),90,'facecolor',[0.3 0.8 0.1])
    bar(200,spont(i,2),90,'facecolor','b')
    bar(300,spont(i,3),90,'facecolor','r')
    

    text(180,32,['ND',int2str(9-i)],opt_title);

    if i==1||i==2
        m=ylabel(sprintf('Spontaneous\nFR,Hz'));
        set(m,opt_label)
    end
    
    set(gca,'xtick',100:100:300,'xticklabel',{'trans','gap','sust'},opt_axis,'fontweight','bold')

    axis([0 400 0 30])

end

saveas(gcf,'S:\user\alexandra\thesis\figures\R4.emf')

%% Fig.R6
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R3_R4_R5_R6.mat','stableCells')
figure
set(gcf,'position',[2357         706         289         217])
subplot('position',[0.15 0.15 0.84 0.84])
bar(stableCells)

axis([0 4 0 25])
set(gca,'xtick',1:3,'xticklabel',{'trans','gap','sust'},'ytick',5:5:20,opt_axis,'fontweight','bold')
m=ylabel('% of ON cells');
set(m,opt_label);

saveas(gcf,'S:\user\alexandra\thesis\figures\R6.emf')


%% Fig.R5
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R3_R4_R5_R6.mat','data','bardata')

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

%% Fig.R12
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R12.mat')

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,1,i),data(:,2,i),'-*')
    hold on
    plot(data(:,3,i),data(:,4,i),'-*r')
    
    axis([-1 1 0 50])
    set(gca,'xtick',[-0.5 0 0.5],opt_axis,'fontweight','bold')
    
    text(-0.1,54,['ND',int2str(9-i)],opt_title);

    if i<4
        m=ylabel('% of ON cells');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('correlation');
        set(m,opt_label)
    end
    line([0 2000],[0 2000],'color','k')

end

subplot(3,3,1,'position',[0.0600    0.7200    0.2870    0.2500])
h=plot(1,1,'r');
hold on
h1=plot(1,1,'b');
m=legend({'ON response','OFF response'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\R12.emf')


%% Fig.R11
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R11.mat','data','bardata')

figure
set(gcf,'position',[1784         201        1101         754])
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])
    
    hold on
    plot(data(:,1,i),'color',[0.3 0.8 0.1],'linewidth',2)
    plot(data(:,2,i),'linewidth',2)
    plot(data(:,3,i),'r','linewidth',2)
    bar(4650,bardata(1,i),180,'facecolor',[0.3 0.8 0.1])
    bar(4850,bardata(2,i),180,'facecolor','b')
    bar(5050,bardata(3,i),180,'facecolor','r')

    axis([0 5250 -20 55])
    
    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500],'xticklabel',{'0.5','2.5'},opt_axis)
    if i==1||i==2
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
    line([1 500],[50,50],'color','k')
    line([500 2500],[55,55]-10,'color','k')
    line([2500 4500],[50,50],'color','k')
    line([500 500],[55,50]-5,'color','k')
    line([2500 2500],[55,50]-5,'color','k')
end
subplot(3,3,1,'position',sizes_1)
h=plot(1,1,'color',[0.3 0.8 0.1]);
hold on
h1=plot(1,1,'b');
h2=plot(1,1,'r');
m=legend({'transient','gap-like','sustained'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
delete(h2)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\R11.emf')


%% Fig.R11a
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R11.mat')

figure
set(gcf,'position',[2162         522         606         293])
subplot('position',[0.08 0.10 0.91 0.89])
a=bar(bardata','hist');
set(a(1),'facecolor',[0.3 0.8 0.1])
set(a(2),'facecolor','b')
set(a(3),'facecolor','r')
set(gca,'xtick',1:8,'xticklabel',{'ND8','ND7','ND6','ND5','ND4','ND3','ND2','ND1'},'ytick',10:20:70,opt_axis)
m=ylabel('%');
set(m,opt_label);
axis([0.3 8.7 0 75])
m=legend({'trans','gap','sust'},'location','northwest');
set(m,opt_leg)
legend boxoff

saveas(gcf,'S:\user\alexandra\thesis\figures\R11a.emf')


%% Fig.R8
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R8.mat')

figure
set(gcf,'position',[1784         201        1101         754])
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])
    
    plot(data{i}(1,:),data{i}(2,:),'o')
    axis([-50 50 -0 2])
    
    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[-25,0,25],'xticklabel',{'-25','0','25'},opt_axis)
    if i==1||i==2
        m=ylabel('amplitude ratio');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('latency delta, ms');
        set(m,opt_label)
    end

    line([0,0],[0,2],'color','k')
    line([-50 50],[1,1],'color','k')
end

subplot(3,3,1,'position',sizes_1)
h1=plot(1,1,'ob');
m=legend({'ON cells'});
set(m,opt_title)
legend('boxoff')
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\R8.emf')

%% Fig.R9
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R9.mat')

figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,i,1),data(:,i,2),'o','markersize',2)
    
    axis([0 2000 0 2000])
    
    text(850,2200,['ND',int2str(9-i)],opt_title);
%     m=title(['ND',int2str(9-i)]);
%     set(m,opt_title);
    set(gca,'xtick',[500,1500],'ytick',[500,1500],opt_axis)
    if i<4
        m=ylabel('Transiency White');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('Transiency Black');
        set(m,opt_label)
    end
    line([0 2000],[0 2000],'color','k')
    p(i)=ranksum(data(:,i,1),data(:,i,2));
end

saveas(gcf,'S:\user\alexandra\thesis\figures\R9.emf')


%% Fig.R10
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_R10.mat')
figure
set(gcf,'position',[ 2210         386         650         436])
for i=1:8
    subplot('position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,i,2),data(:,i,1),'o','markersize',2)
    
    axis([0 2000 0 2000])
    
    text(850,2200,['ND',int2str(9-i)],opt_title);
%     m=title(['ND',int2str(9-i)]);
%     set(m,opt_title);
    set(gca,'xtick',[500,1500],'ytick',[500,1500],opt_axis)
    if i<4
        m=ylabel('Transiency White');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('Transiency Black');
        set(m,opt_label)
    end
    line([0 2000],[0 2000],'color','k')
    p(i)=ranksum(data(:,i,1),data(:,i,2));
end

saveas(gcf,'S:\user\alexandra\thesis\figures\R10.emf')


%% Fig.R13

clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_Z6_alt.mat')

for i=1:8
    % cut out noise
%     cc=abs(data{i}(1,:))<55&data{i}(2,:)<2&data{i}(2,:)>0.2;
    cc=abs(data{i}(1,:))<25&data{i}(2,:)<1.5&data{i}(2,:)>0.5;
    
    % amplitudes ratio
    tmp=data{i}(2,cc);
    k(i)=nanmean(tmp);
    ks(i)=nanstd(tmp)/sqrt(size(tmp,2));
    amp_black{i}=tmp;
    [n(i),b(i)]=ttest(tmp-1);
    
    r(i)=kstest(tmp);
    
    % latency difference
    tmp1=data{i}(1,cc);
    t(i)=nanmean(tmp1);
    ts(i)=nanstd(tmp1)/sqrt(size(tmp1,2));
    lat_black{i}=tmp1;
    [m(i),a(i)]=ttest(tmp1);  

    r1(i)=kstest(tmp1);
    
    % correlation of latency and amplitude
    kp(i)=corr(tmp',tmp1');
end

figure
set(gcf,'position',[1973         599         866         273])

sb1=subplot('position',[0.06 0.15 0.43 0.83])
errorbar(k,ks,'b')
hold on

sb2=subplot('position',[0.56 0.15 0.43 0.83])
errorbar(t,ts,'b')
hold on

load('S:\data\alexandra\MEA_data\analysis\data_for_R8.mat')

for i=1:8
    % cut out noise
%     cc=abs(data{i}(1,:))<55&data{i}(2,:)<2&data{i}(2,:)>0.2;
    cc=abs(data{i}(1,:))<25&data{i}(2,:)<1.5&data{i}(2,:)>0.5;
    
    % amplitudes ratio
    tmp=data{i}(2,cc);
    k(i)=nanmean(tmp);
    ks(i)=nanstd(tmp)/sqrt(size(tmp,2));
    amp_white{i}=tmp;
    [n(i),b(i)]=ttest(tmp-1);
    
    r(i)=kstest(tmp);
    
    % latency difference
    tmp1=data{i}(1,cc);
    t(i)=nanmean(tmp1);
    ts(i)=nanstd(tmp1)/sqrt(size(tmp1,2));
    lat_white{i}=tmp1;
    [m(i),a(i)]=ttest(tmp1);  
    r1(i)=kstest(tmp1);
    
    % correlation of latency and amplitude
    kp(i)=corr(tmp',tmp1');
end

for i=1:8
    amp(i)=ranksum(amp_black{i},amp_white{i});
    lat(i)=ranksum(lat_black{i},lat_white{i});
end


subplot(sb1)
errorbar(k,ks,'r')
m=legend({'OFF cells OFF response','ON cells ON response'});
set(m,opt_leg);
legend boxoff
for i=1:8
    if amp(i)<0.001
        plot(i-0.1,k(i)+ks(i)+0.07,'*k')
        plot(i+0.1,k(i)+ks(i)+0.07,'*k')
    elseif  amp(i)<0.05
        plot(i,k(i)+ks(i)+0.7,'*k')
    end
end
line([0 9],[1,1],'color','b')
axis([0 9 0.8 1.4])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Amplitude ratio');
set(m,opt_label);

subplot(sb2)
errorbar(t,ts,'r')
m=legend({'OFF cells OFF response','ON cells ON response'});
set(m,opt_leg);
legend boxoff
for i=1:8
    if lat(i)<0.001
        plot(i-0.1,t(i)+ts(i)+0.7,'*k')
        plot(i+0.1,t(i)+ts(i)+0.7,'*k')
    elseif  lat(i)<0.05
        plot(i,t(i)+ts(i)+0.7,'*k')     
    end
end
line([0 9],[0,0],'color','k')
axis([0 9 -15 5])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'},opt_axis)
m=xlabel('ND');
set(m,opt_label);
m=ylabel('Latency difference,ms');
set(m,opt_label);


saveas(gcf,'S:\user\alexandra\thesis\figures\R13.emf')

%% Fig.P
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_P.mat')

figure
set(gcf,'position',[1784         201        1101         754])
for i=1:8
    subplot(3,3,plnmr(i),'position',[corn_posx(i) corn_posy(i) size_x size_y])

    plot(data(:,1,i),'r');
    hold on
    plot(data(:,2,i),'b');
    
    plot(4701:9200,data(:,3,i),'r');
    plot(4701:9200,data(:,4,i),'b');

    m=title(['ND',int2str(9-i)]);
    set(m,opt_title);
    set(gca,'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2.5','0.5','2.5'},opt_axis)
    if i==1||i==2
        m=ylabel('Firing rate,Hz');
        set(m,opt_label)
    end
    if i==2||i==5||i==8
        m=xlabel('time,s');
        set(m,opt_label)
    end
    
    axis([1 9200 -20 65])
    
    line([1 500],[55,55],'color','k')
    line([500 2500],[60,60],'color','k')
    line([2500 4500],[55,55],'color','k')
    line([500 500],[55,60],'color','k')
    line([2500 2500],[55,60],'color','k')
    
    line([1 500]+4700,[55,55],'color','k')
    line([500 2500]+4700,[50,50],'color','k')
    line([2500 4500]+4700,[55,55],'color','k')
    line([500 500]+4700,[50,55],'color','k')
    line([2500 2500]+4700,[50,55],'color','k')
end
subplot(3,3,1,'position',sizes_1)
h=plot(1,1,'r');
hold on
h1=plot(1,1,'b');
m=legend({'ON cells','OFF cells'});
set(m,opt_title)
legend('boxoff')
delete(h)
delete(h1)
set(gca,'Visible','off')

saveas(gcf,'S:\user\alexandra\thesis\figures\P.emf')


%% Fig.O
clear
load('S:\data\alexandra\MEA_data\analysis\plot_opts')
load('S:\data\alexandra\MEA_data\analysis\data_for_O.mat')

figure
set(gcf,'position',[2170         351         593         364])
set(gcf,'color',[1 1 1])
line_coord(1)=3950;
line_coord(2)=3150;
line_coord(3)=5450;
line_coord(4)=4850;
line_coord(5)=13050;

for i=1:3
    p=subplot('Position',[0.38 0.08+0.298*(i-1)+0.11 0.58 0.16]);
    hold on
    line([line_coord(1) line_coord(1)],[0 36],'color',[1 0.8 0.8]*0.9,'linewidth',4)
    line([line_coord(2) line_coord(2)],[0 36],'color',[0.8 1 0.8]*0.9,'linewidth',5)
    line([line_coord(3) line_coord(3)],[0 36],'color',[0.8 1 0.8]*0.9,'linewidth',7)
    line([line_coord(4) line_coord(4)],[0 36],'color',[1 0.8 0.8]*0.9,'linewidth',3)
    line([line_coord(5) line_coord(5)],[0 36],'color',[1 0.8 0.8]*0.9,'linewidth',5)
    rasterplot(dataRaster{i},dataRaster_supp(i),26000,p);


    set(gca,'xtick',0,'xticklabel','')
    set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])
    axis([0 15675 0 27])
    
    p=subplot('Position',[0.38 0.08+0.298*(i-1) 0.58 0.11]);
    plot(conv_acc(i,1:25000),'k','linewidth',2)
%     axis tight
    set(gca,'xtick',0,'xticklabel','')
    set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])
    axis([0 15675 0 130])
    set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])
    
    if i==1
        subplot('position',[0.38 0.01 0.074 0.07])
        line([0 2],[0,0],'color','k','linewidth',2)
        text(0.5,0.3,'2s',opt_text)
        axis([0 2 -0.1 1])
        set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])
    end
    
end


subplot('Position',[0.01 0.02 0.34 0.97])
hold on
beg1=1;
beg2=4701;
hline(1)=line([beg1+500 beg1+500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
hline(2)=line([beg1+2500 beg1+2500],[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
hline(3)=line([beg1+500 beg1+500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);
hline(4)=line([beg1+2500 beg1+2500]+beg2-beg1,[-20,1000],'color',[1 1 1]*0.5,'linewidth',2);

const=0;
stim_height=5;
for i=1:3

    line([0 beg2+4500],[const const],'color',[1 1 1]*0.5,'linewidth',2);    
    plot(beg1:beg1+4499,data(:,1,i)+const,'k','linewidth',2);   
    plot(beg2:beg2+4499,data(:,2,i)+const,'k','linewidth',2);    
    text(beg2+3400,const+50,['ND',int2str(7-i)],opt_text,'fontsize',14)
    const=const+100;    
    
end
const=const+5;
axis([1 beg2+4500 -20 const])

k=-15;
for i=1:4
    set(hline(i),'ydata',[k 300])
end

line([beg1 beg1+500],[0,0]+k,'color','k','linewidth',2);
line([beg1+500 beg1+2500],[stim_height,stim_height]+k,'color','k','linewidth',2);
line([beg1+2500 beg1+4500],[0,0]+k,'color','k','linewidth',2);
line([beg1+500 beg1+500],[-1,stim_height+1]+k,'color','k','linewidth',2);
line([beg1+2500 beg1+2500],[-1,stim_height+1]+k,'color','k','linewidth',2);

line([beg1 beg1+500]+beg2-beg1,[0,0]+k,'color','k','linewidth',2);
line([beg1+500 beg1+2500]+beg2-beg1,[-stim_height,-stim_height]+k,'color','k','linewidth',2);
line([beg1+2500 beg1+4500]+beg2-beg1,[0,0]+k,'color','k','linewidth',2);
line([beg1+500 beg1+500]+beg2-beg1,[-1.5-stim_height,1.5]+k,'color','k','linewidth',2);
line([beg1+2500 beg1+2500]+beg2-beg1,[-1.5-stim_height,1.5]+k,'color','k','linewidth',2);

line([3200 3200],[120,160],'color','k','linewidth',2);
text(3350,140,'40Hz',opt_text)

set(gcf,'color',[1 1 1])
set(gca,'ycolor',[1 1 1],'xcolor',[1 1 1])


saveas(gcf,'S:\user\alexandra\thesis\figures\O.emf')

