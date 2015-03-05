%% Fig.K0
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
figure

data=zeros(4,8,6);
clear m v
for i=1:8
    m(i)=mean(frMean_spont(onOff>0,i));
    m1(i)=std(frMean_spont(onOff>0,i))/sqrt(sum(onOff>0));
    v(i)=mean(frMean_spont(onOff<0,i));
    v1(i)=std(frMean_spont(onOff<0,i))/sqrt(sum(onOff<0));
end
data(1,:,1)=m;
data(2,:,1)=m1;
data(3,:,1)=v;
data(4,:,1)=v1;

subplot(2,3,1)
errorbar(m,m1,'r')
hold on
errorbar(v,v1,'b')
title('Mean Spont FR')
axis([0 9 0 30])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('Hz')
legend('ON','OFF')

clear m v
for i=1:8
    m(i)=mean(frSTD_spont(onOff>0,i));
    m1(i)=std(frSTD_spont(onOff>0,i))/sqrt(sum(onOff>0));
    v(i)=mean(frSTD_spont(onOff<0,i));
    v1(i)=std(frSTD_spont(onOff<0,i))/sqrt(sum(onOff<0));
end
data(1,:,4)=m;
data(2,:,4)=m1;
data(3,:,4)=v;
data(4,:,4)=v1;

subplot(2,3,4)
errorbar(m,m1,'r')
hold on
errorbar(v,v1,'b')
title('Variance Spont FR')
axis([0 9 0 20])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('Hz')

clear m v
for i=1:8
    m(i)=mean(frMean_HC(onOff>0,i));
    m1(i)=std(frMean_HC(onOff>0,i))/sqrt(sum(onOff>0));
    v(i)=mean(frMean_HC(onOff<0,i));
    v1(i)=std(frMean_HC(onOff<0,i))/sqrt(sum(onOff<0));
end
data(1,:,2)=m;
data(2,:,2)=m1;
data(3,:,2)=v;
data(4,:,2)=v1;

subplot(2,3,2)
errorbar(m,m1,'r')
hold on
errorbar(v,v1,'b')
title('Mean HC FR')
axis([0 9 0 30])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('Hz')

clear m v
for i=1:8
    m(i)=mean(frSTD_HC(onOff>0,i));
    m1(i)=std(frSTD_HC(onOff>0,i))/sqrt(sum(onOff>0));
    v(i)=mean(frSTD_HC(onOff<0,i));
    v1(i)=std(frSTD_HC(onOff<0,i))/sqrt(sum(onOff<0));
end
data(1,:,5)=m;
data(2,:,5)=m1;
data(3,:,5)=v;
data(4,:,5)=v1;

subplot(2,3,5)
errorbar(m,m1,'r')
hold on
errorbar(v,v1,'b')
title('Variance HC FR')
axis([0 9 0 20])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('Hz')


clear m v
for i=1:8
    m(i)=mean(frMean_LC(onOff>0&bothContr,i));
    m1(i)=std(frMean_LC(onOff>0&bothContr,i))/sqrt(sum(onOff>0&bothContr));
    v(i)=mean(frMean_LC(onOff<0&bothContr,i));
    v1(i)=std(frMean_LC(onOff<0&bothContr,i))/sqrt(sum(onOff<0&bothContr));
end
data(1,:,3)=m;
data(2,:,3)=m1;
data(3,:,3)=v;
data(4,:,3)=v1;

subplot(2,3,3)
errorbar(m,m1,'r')
hold on
errorbar(v,v1,'b')
title('Mean LC FR')
axis([0 9 0 30])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('Hz')


clear m v
for i=1:8
    m(i)=mean(frSTD_LC(onOff>0&bothContr,i));
    m1(i)=std(frSTD_LC(onOff>0&bothContr,i))/sqrt(sum(onOff>0&bothContr));
    v(i)=mean(frSTD_LC(onOff<0&bothContr,i));
    v1(i)=std(frSTD_LC(onOff<0&bothContr,i))/sqrt(sum(onOff<0&bothContr));
end
data(1,:,6)=m;
data(2,:,6)=m1;
data(3,:,6)=v;
data(4,:,6)=v1;
subplot(2,3,6)
errorbar(m,m1,'r')
hold on
errorbar(v,v1,'b')
title('Variance LC FR')
axis([0 9 0 20])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('Hz')

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_K0','data')



%% Fig.K01
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','frMean_spont','onOff','frSTD_spont','frMean_HC',...
    'frSTD_HC','frMean_LC','frSTD_LC','bothContr')
figure
set(gcf,'position',[1743         443         973         492])
clear m v
for i=1:8
    s(i)=mean(frMean_spont(onOff>0,i));
    s1(i)=std(frMean_spont(onOff>0,i))/sqrt(sum(onOff>0));
    m(i)=mean(frMean_HC(onOff>0,i));
    m1(i)=std(frMean_HC(onOff>0,i))/sqrt(sum(onOff>0));
    v(i)=mean(frMean_LC(onOff>0&bothContr,i));
    v1(i)=std(frMean_LC(onOff>0&bothContr,i))/sqrt(sum(onOff>0&bothContr));
end
subplot(1,2,1)
errorbar(m,m1,'m')
hold on
errorbar(v,v1,'c')
errorbar(s,s1,'k')
title('Mean firing rate, ON cells')
axis([0 9 0 30])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND','fontweight','bold','fontsize',12)
ylabel('Hz','fontweight','bold','fontsize',12)
legend({'High contrast','Low contrast','spontaneous'},'fontsize',12)


clear m v
for i=1:8
    s(i)=mean(frMean_spont(onOff<0,i));
    s1(i)=std(frMean_spont(onOff<0,i))/sqrt(sum(onOff<0));
    m(i)=mean(frMean_HC(onOff<0,i));
    m1(i)=std(frMean_HC(onOff<0,i))/sqrt(sum(onOff<0));
    v(i)=mean(frMean_LC(onOff<0&bothContr,i));
    v1(i)=std(frMean_LC(onOff<0&bothContr,i))/sqrt(sum(onOff<0&bothContr));
end
subplot(1,2,2)
errorbar(m,m1,'m')
hold on
errorbar(v,v1,'c')
errorbar(s,s1,'k')
title('Mean firing rate, OFF cells')
axis([0 9 0 30])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND','fontweight','bold','fontsize',12)
ylabel('Hz','fontweight','bold','fontsize',12)
legend({'High contrast','Low contrast','spontaneous'},'fontsize',12)
