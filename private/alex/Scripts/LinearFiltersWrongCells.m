
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')

figure
nds='87654321'
ons=onOff>0&bothContr>0;
for i=1:8
    subplot(2,4,i)
    k=((frMean_HC(:,i)-frMean_spont(:,i))<-0.2*mean(frMean_spont(:,i)))';
    t=reshape(LC(:,i,k&ons),500,sum(k&ons));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2))
    hold on
    t=reshape(HC(:,i,k&ons),500,sum(k&ons));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2),'r')
    legend({'ON wrong low','ON wrong high'})
    line([0,500],[0,0],'color','k')
    title([nds(i),'  ',int2str(sum(k&ons))])
end

figure
for i=1:8
    subplot(2,4,i)
    k=((frMean_HC(:,i)-frMean_spont(:,i))>0.2*mean(frMean_spont(:,i)))';
    t=reshape(LC(:,i,k&ons),500,sum(k&ons));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2))
    hold on
    t=reshape(HC(:,i,k&ons),500,sum(k&ons));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2),'r')
    legend({'ON reg low','ON reg high'})
    line([0,500],[0,0],'color','k')
    title([nds(i),'  ',int2str(sum(k&ons))])
end

figure
ons=onOff<0&bothContr>0;
for i=1:8
    subplot(2,4,i)
    k=((frMean_HC(:,i)-frMean_spont(:,i))>0.2*mean(frMean_spont(:,i)))';
    t=reshape(LC(:,i,k&ons),500,sum(k&ons));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2))
    hold on
    t=reshape(HC(:,i,k&ons),500,sum(k&ons));
%     t=t-repmat(mean(t),size(t,1),1);
    t=t./repmat(sum(abs(t)),size(t,1),1);
%     plot(t,'b')
    plot(nanmean(t,2),'r')
    legend({'OFF low','OFF high'})
    line([0,500],[0,0],'color','k')
    title([nds(i),'  ',int2str(sum(k&ons))])
end

%% combined

figure
nds='87654321'
ons=onOff>0&bothContr>0;
offs=onOff<0&bothContr>0;
for i=1:8
    
    subplot(2,4,i)
    
    % ON wrong
    k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))<-0.05)';
    t=reshape(LC(:,i,k&ons),500,sum(k&ons));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2))
    hold on
    t=reshape(HC(:,i,k&ons),500,sum(k&ons));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'r')
    
    
    % ON reg
    k1=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))>0.05)';
    t=reshape(LC(:,i,k1&ons),500,sum(k1&ons));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'c')
    t=reshape(HC(:,i,k1&ons),500,sum(k1&ons));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'m')
    
    
    % OFF reg
    
    k2=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))>0.05)';
    t=reshape(LC(:,i,k2&offs),500,sum(k2&offs));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'k')
    t=reshape(HC(:,i,k2&offs),500,sum(k2&offs));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'color',[0.3 0.8 0.1])
    
    % OFF very reg
    k3=((frMean_HC(:,i)-frMean_spont(:,i))>=0.1*mean(frMean_spont(:,i)))';
    t=reshape(LC(:,i,k3&offs),500,sum(k3&offs));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(-nanmean(t,2),'k','linewidth',2)
    t=reshape(HC(:,i,k3&offs),500,sum(k3&offs));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(-nanmean(t,2),'color',[0.3 0.8 0.1],'linewidth',2)
    
    
    legend({['ON wrong low ', int2str(sum(k&ons))],['ON wrong high ', int2str(sum(k&ons))],...
        ['ON reg low ', int2str(sum(k1&ons))],['ON reg high ', int2str(sum(k1&ons))],...
        ['OFF low ', int2str(sum(k2&offs))],['OFF high ', int2str(sum(k2&offs))],...
        ['OFF very low ', int2str(sum(k3&offs))], ['OFF very high ', int2str(sum(k3&offs))]})
    line([0,500],[0,0],'color','k')
    title([nds(i),'  ',int2str(sum(k&ons))])
end

figure

for i=1:8
    subplot(2,4,i)
    k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))<-0.05)';
    a=(frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i));
    b=(zc_HC(i,:)-zc_LC(i,:))';
    mm(i,1)=mean(b(onOff>0&k&abs(b')<20&bothContr>0));
    mm(i,2)=mean(b(onOff>0&~k&abs(b')<20&bothContr>0));
    mm(i,3)=mean(b(onOff<0&abs(b')<20&bothContr>0));
    plot(a(onOff>0&k),b(onOff>0&k),'.g')
    hold on
    plot(a(onOff>0&~k&bothContr>0),b(onOff>0&~k&bothContr>0),'.r')
    plot(a(onOff<0&bothContr>0),b(onOff<0&bothContr>0),'.b')
    axis([-1 1 -30 30])
    line([-3 3],[0 0],'color','k')
    line([0 0],[-50 50],'color','k')
end





%% chirp


figure
nds='87654321'
ons=onOff>0;
offs=onOff<0;
for i=1:8
    
    subplot(4,2,i)
    
    % ON wrong
   k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))<-0.05)';
    
    a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'&ons),2);
    a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);    
    a=mean(a,2);
    plot(a,'r')    
    hold on   
    
    for indd=1:length(maxpos)
        [maxx(indd,i,1) indmax(indd,i,1)]=max(a(maxpos(indd)-100:maxpos(indd)+100));
        indmax(indd,i,1)=indmax(indd,i,1)+maxpos(indd)-100-1;
        [minn(indd,i,1) indmin(indd,i,1)]=min(a(minneg(indd)-100:minneg(indd)+100));
        indmin(indd,i,1)=indmin(indd,i,1)+minneg(indd)-100-1;
    end
    
    
    
    % ON reg
    k1=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))>0.05)';
    a=mean(chirp(:,(i-1)*4+1:i*4,k1&real_chirp'&ons),2);
    a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);    
    a=mean(a,2);
    plot(a,'b') 
    
    for indd=1:length(maxpos)
        [maxx(indd,i,2) indmax(indd,i,2)]=max(a(maxpos(indd)-100:maxpos(indd)+100));
        indmax(indd,i,2)=indmax(indd,i,2)+maxpos(indd)-100-1;
        [minn(indd,i,2) indmin(indd,i,2)]=min(a(minneg(indd)-100:minneg(indd)+100));
        indmin(indd,i,2)=indmin(indd,i,2)+minneg(indd)-100-1;
    end
    
    
    
    
    % OFF reg
    
    a=mean(chirp(:,(i-1)*4+1:i*4,real_chirp'&offs),2);
    a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);    
    a=mean(a,2);
    plot(a,'color',[0.3 0.8 0.1])
    for indd=1:length(maxpos)
        [minn(indd,i,3) indmin(indd,i,3)]=min(a(maxpos(indd)-200:maxpos(indd)+200));
        indmin(indd,i,3)=indmin(indd,i,3)+maxpos(indd)-200-1;
        [maxx(indd,i,3) indmax(indd,i,3)]=max(a(minneg(indd)-200:minneg(indd)+200));
        indmax(indd,i,3)=indmax(indd,i,3)+minneg(indd)-200-1;
    end
   
    
    legend({['ON wrong ', int2str(sum(k&ons&real_chirp'))],...
        ['ON reg ', int2str(sum(k1&ons&real_chirp'))],...
        ['OFF ', int2str(sum(real_chirp'&offs))]})
    line([0,500],[0,0],'color','k')
    title([nds(i)])
end



%contrast response function
figure
for i=1:8
    subplot(4,2,i)
    plot([flicker(minneg(end:-1:1))-30; flicker(maxpos)-30],[minn(end:-1:1,i,1); maxx(:,i,1)],'r')
    hold on
    plot([flicker(minneg(end:-1:1))-30; flicker(maxpos)-30],[minn(end:-1:1,i,2); maxx(:,i,2)])
    plot([flicker(minneg(end:-1:1))-30; flicker(maxpos)-30],[maxx(end:-1:1,i,3); minn(:,i,3)],'g')
    legend({'wrong cells','reg cells','off'},'location','southeast')
    line([-30 30],[0,0],'color','k')
    line([0,0],[-30 30],'color','k')
    title(nds(i))
%     axis([-10 10 0 40])
end



figure
nds='87654321'
ons=onOff>0;
offs=onOff<0;
i=7

% ON wrong
k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))<-0.05)';

a=mean(chirp(:,(i-1)*4+1:i*4,k&real_chirp'&ons),2);
a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);
a=mean(a,2);
plot(a,'r')
hold on

% ON reg
k1=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))>0.05)';
a=mean(chirp(:,(i-1)*4+1:i*4,k1&real_chirp'&ons),2);
a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);
a=mean(a,2);
plot(a,'b')

% OFF reg

a=mean(chirp(:,(i-1)*4+1:i*4,real_chirp'&offs),2);
a=reshape(a,23000,size(a,3));
%     a=a-repmat(mean(a(500:2500,:)),23000,1);
a=mean(a,2);
plot(a,'color',[0.3 0.8 0.1])


legend({['ON wrong ', int2str(sum(k&ons&real_chirp'))],...
    ['ON reg ', int2str(sum(k1&ons&real_chirp'))],...
    ['OFF ', int2str(sum(real_chirp'&offs))]})
line([0,500],[0,0],'color','k')
title([nds(i)])



%% quick

figure
nds='87654321'
ons=onOff>0;
for i=1:8
    
    subplot(2,4,i)
    
    % ON wrong
    k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))<-0.05)';
    t=wf{i}(:,k&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(nanmean(t,2),'r')
    hold on
    t=bf{i}(:,k&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(4701:9200,nanmean(t,2),'r')
    
    % ON reg
    k1=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))>0.05)';
    t=wf{i}(:,k1&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(nanmean(t,2))
    t=bf{i}(:,k1&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(4701:9200,nanmean(t,2))
    
    
    legend({['ON wrong ', int2str(sum(k&ons))],['ON reg ', int2str(sum(k1&ons))]})
    line([0,500],[0,0],'color','k')
    title([nds(i)])
end






figure
nds='87654321'

for i=1:8
    
    subplot(2,4,i)
    
    % ON wrong
    k=(frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i))./(frMean_HC(onOff>0,i)+frMean_spont(onOff>0,i));
    plot(sort(k),'.r')
    hold on
    k=(frMean_HC(onOff<0,i)-frMean_spont(onOff<0,i))./(frMean_HC(onOff<0,i)+frMean_spont(onOff<0,i));
    plot(sort(k),'.')
    line([0,250],[0,0],'color','k')
    title([nds(i)])
end







figure
nds='87654321'
ons=onOff>0;
for i=1:8
    
    subplot(2,4,i)
    
    % ON wrong
    k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))<-0.05)';
    t=frMean_spont(k&onOff>0,i);
    plot(t,'.r')
    hold on
    
    % ON reg
    k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))>0.05)';
    t=frMean_spont(k&onOff>0,i);
    plot(t,'.')
    title([nds(i)])
end





figure
nds='87654321'
ons=onOff>0;
for i=1:8
    
    subplot(2,4,i)
    
    % ON wrong
    k=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))<-0.05)';
    t=wf{i}(:,k&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(nanmean(t,2),'r')
    hold on
    t=bf{i}(:,k&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(4701:9200,nanmean(t,2),'r')
    
    % ON reg
    k1=((frMean_HC(:,i)-frMean_spont(:,i))./(frMean_HC(:,i)+frMean_spont(:,i))>0.05)';
    t=wf{i}(:,k1&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(nanmean(t,2))
    t=bf{i}(:,k1&ons);
    t=t-repmat(mean(t(50:450,:)),4500,1);
    plot(4701:9200,nanmean(t,2))
    
    
    legend({['ON wrong ', int2str(sum(k&ons))],['ON reg ', int2str(sum(k1&ons))]})
    line([0,500],[0,0],'color','k')
    title([nds(i)])
end



%% OFF cells

figure
cc=1;tt=1;
for code=code_list([1 5 2 6 3 4])
    cc=tt;
    for i=3:5
        exx=int2str(length(unique(exp_codes(a==code&bothContr'))));
        subplot(3,6,cc)
        t=HC(:,i,a==code&bothContr');   
        t=t./repmat(sum(abs(t)),size(t,1),1);
        plot(nanmean(t,3)*1000,'color','b','linewidth',2)
        hold on
        t=LC(:,i,a==code&bothContr');
        t=t./repmat(sum(abs(t)),size(t,1),1);
        plot(nanmean(t,3)*1000,'color','r','linewidth',2)        
        
        title(['ND',nds(i),' n=',int2str(sum(a==code)), ' in ',exx,' experiments'])
        axis([0 500 -10 10])
        cc=cc+6;
    end
    tt=tt+1;
end




figure
cc=1;tt=1;
col='rgbmck'
for code=code_list([1 5 2 6 3 4])
    cc=1;
    for i=3:5
        subplot(1,3,cc)
        t=HC(:,i,a==code);   
        t=t./repmat(sum(abs(t)),size(t,1),1);
        plot(nanmean(t,3)*1000,'color',col(tt),'linewidth',2)
        hold on     
        
        title(['ND',nds(i)])
        axis([0 500 -12 12])
        cc=cc+1;
    end
    tt=tt+1;
end


%% ON cells

% Fig R6
figure
for i=1:8
    subplot(2,4,i)
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    
    ff=mean(tmp(800:2450,:))';
    ff=ff<1;
    kf=mean(tmp(500:700,:))';
    fi=mean(tmp(680:820,:))';
    ti=mean(tmp(880:1000,:))';
    si=mean(tmp(2000:2450,:))';
    si=si>1;
    ti=ti>fi;
    fi=fi<(kf-3);
    kf=kf>3;
    trans=kf&ff;
    inh=fi&kf&(~ff)&ti;
    sus=si&kf&(~ff)&(~inh);
    m0=find(trans);
    m=find(inh);
    m1=find(sus);
    m2=ones(size(kf));
    m2([m0; m; m1])=0;
    m2=find(m2);
    
    hold on
    t=HC(:,i,onOff>0);
    t=t./repmat(sum(abs(t)),size(t,1),1);
    
    plot(nanmean(t(:,1,m0),3)*1000,'color',[0.3 0.8 0.1],'linewidth',2)
    plot(nanmean(t(:,1,m),3)*1000,'linewidth',2)
    plot(nanmean(t(:,1,m1),3)*1000,'color','r','linewidth',2)
    plot(nanmean(t(:,1,m2),3)*1000,'color','k','linewidth',2)
    
    axis([0 500 -5 15])
    title(['ND',nds(i),' trans ',int2str(round(length(m0)/2.49)),'% inh ',int2str(round(length(m)/2.49)),'% sust ',int2str(round(length(m1)/2.49)),'%'])
    legend('trans','inh','sust')

end




% Fig R6
figure
for i=1:8
    subplot(2,4,i)
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    
    ff=mean(tmp(800:2450,:))';
    ff=ff<1;
    kf=mean(tmp(500:700,:))';
    fi=mean(tmp(680:820,:))';
    ti=mean(tmp(880:1000,:))';
    si=mean(tmp(2000:2450,:))';
    si=si>1;
    ti=ti>fi;
    fi=fi<(kf-3);
    kf=kf>3;
    trans=kf&ff;
    inh=fi&kf&(~ff)&ti;
    sus=si&kf&(~ff)&(~inh);
    m0=find(trans);
    m=find(inh);
    m1=find(sus);
    m2=ones(size(kf));
    m2([m0; m; m1])=0;
    m2=find(m2);
    
    k=((frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i))./(frMean_HC(onOff>0,i)+frMean_spont(onOff>0,i))<-0.05)';
    a=[sum(k(m0)) sum(k(m)) sum(k(m1)) sum(k(m2))];
    
    bar(a)
    
    title(['ND',nds(i),' trans ',num2str(a(1)/sum(k)),'% inh ',num2str(a(2)/sum(k)),'% sust ',num2str(a(3)/sum(k)),'% rest ',num2str(a(4)/sum(k)),'%'])
end



%% ZC

tmp=zc_HC-zc_LC;
cc=1;
for i=1:8
    subplot(4,4,cc)
    cond=zc_LC(i,:)>0&bothContr;
    hist([-200,100, zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0)],200);
    line([0,0],[0,20],'color','r')
    axis([-40 20 0 25])
    subplot(4,4,cc+1)
    hist([-200,100, zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0)],200);
    axis([-40 20 0 25])
    line([0,0],[0,20],'color','r')
    cc=cc+2;
end

figure
tmp=peak_HC-peak_LC;
cc=1;
for i=1:8
    subplot(4,4,cc)
    cond=peak_LC(i,:)>0&bothContr;
    hist([-200,100, peak_HC(i,cond&onOff>0)-peak_LC(i,cond&onOff>0)],200);
    line([0,0],[0,20],'color','r')
    axis([-40 20 0 25])
    subplot(4,4,cc+1)
    hist([-200,100, peak_HC(i,cond&onOff<0)-peak_LC(i,cond&onOff<0)],200);
    axis([-40 20 0 25])
    line([0,0],[0,20],'color','r')
    cc=cc+2;
end





clear a b c d
figure
for i=1:8
    tmp=zc_HC-zc_LC;
    cond=zc_LC(i,:)>0&bothContr&tmp(i,:)<-3&abs(tmp(i,:))<15;
    n=length(zc_HC(i,cond&onOff>0));
    n1=length(zc_HC(i,cond&onOff<0));
    a(i)=mean(zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0));
    b(i)=std(zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0))/sqrt(n);
    c(i)=mean(zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0));
    d(i)=std(zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0))/sqrt(n1);
    ons(i)=n;
    offs(i)=n1;
end
errorbar(a,b,'r')
hold on
errorbar(1.1:1:8.1,c,d,'b')


clear a b c d
figure
for i=1:8
    tmp=zc_HC-zc_LC;
    cond=zc_LC(i,:)>0&bothContr&tmp(i,:)>3&abs(tmp(i,:))<15;
    n=length(zc_HC(i,cond&onOff>0));
    n1=length(zc_HC(i,cond&onOff<0));
    a(i)=mean(zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0));
    b(i)=std(zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0))/sqrt(n);
    c(i)=mean(zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0));
    d(i)=std(zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0))/sqrt(n1);
    ons(i)=n;
    offs(i)=n1;
end
errorbar(a,b,'r')
hold on
errorbar(1.1:1:8.1,c,d,'b')






figure
for i=1:8
    a(i)=mean(zc_HC(i,zc_HC(i,:)>0&onOff>0));
    b(i)=std(zc_HC(i,zc_HC(i,:)>0))/sqrt(sum(zc_HC(i,:)>0&onOff>0));
    c(i)=mean(zc_HC(i,zc_HC(i,:)>0&onOff<0));
    d(i)=std(zc_HC(i,zc_HC(i,:)>0))/sqrt(sum(zc_HC(i,:)>0&onOff<0));
end
errorbar(a,b,'r')
hold on
errorbar(c,d,'b')

for i=1:8
    a(i)=mean(zc_LC(i,zc_LC(i,:)>0&onOff>0&bothContr));
    b(i)=std(zc_LC(i,zc_LC(i,:)>0))/sqrt(sum(zc_LC(i,:)>0&onOff>0&bothContr));
    c(i)=mean(zc_LC(i,zc_LC(i,:)>0&onOff<0&bothContr));
    d(i)=std(zc_LC(i,zc_LC(i,:)>0))/sqrt(sum(zc_LC(i,:)>0&onOff<0&bothContr));
end
errorbar(a,b,'m')
hold on
errorbar(c,d,'c')



% figure
a=peak_HC-peak_LC;
peak_HC(abs(a)>30)=0;
peak_LC(abs(a)>30)=0;
clear a

for i=1:8
    a(i)=mean(peak_HC(i,peak_HC(i,:)>0&onOff>0));
    b(i)=std(peak_HC(i,peak_HC(i,:)>0))/sqrt(sum(peak_HC(i,:)>0&onOff>0));
    c(i)=mean(peak_HC(i,peak_HC(i,:)>0&onOff<0));
    d(i)=std(peak_HC(i,peak_HC(i,:)>0))/sqrt(sum(peak_HC(i,:)>0&onOff<0));
end
errorbar(a,b,'r','linewidth',2)
hold on
errorbar(c,d,'b','linewidth',2)

for i=1:8
    a(i)=mean(peak_LC(i,peak_LC(i,:)>0&onOff>0&bothContr));
    b(i)=std(peak_LC(i,peak_LC(i,:)>0))/sqrt(sum(peak_LC(i,:)>0&onOff>0&bothContr));
    c(i)=mean(peak_LC(i,peak_LC(i,:)>0&onOff<0&bothContr));
    d(i)=std(peak_LC(i,peak_LC(i,:)>0))/sqrt(sum(peak_LC(i,:)>0&onOff<0&bothContr));
end
errorbar(a,b,'m','linewidth',2)
hold on
errorbar(c,d,'c','linewidth',2)

legend({'zc HC ON','zc HC OFF','zc LC ON','zc LC OFF','peak HC ON','peak HC OFF','peak LC ON','peak LC OFF'})



figure
for i=1:8
    subplot(2,4,i)
    hold on
    a=HC(:,i,peak_HC(i,:)>0&onOff>0);
    t=reshape(a,500,size(a,3));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'r','linewidth',2)
    [~,l(i,1)]=max(nanmean(t,2));
    
    a=HC(:,i,peak_HC(i,:)>0&onOff<0);
    t=reshape(a,500,size(a,3));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(-nanmean(t,2),'b','linewidth',2)
    [~,l(i,2)]=max(-nanmean(t,2));
    
    a=LC(:,i,peak_LC(i,:)>0&onOff>0&bothContr);
    t=reshape(a,500,size(a,3));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(nanmean(t,2),'m','linewidth',2)
    [~,l(i,3)]=max(nanmean(t,2));
    
    a=LC(:,i,peak_LC(i,:)>0&onOff<0&bothContr);
    t=reshape(a,500,size(a,3));
    t=t./repmat(sum(abs(t)),size(t,1),1);
    plot(-nanmean(t,2),'c','linewidth',2)
    [~,l(i,4)]=max(-nanmean(t,2));
end

figure
plot(l)
legend({'peak HC ON','peak HC OFF','peak LC ON','peak LC OFF'})
