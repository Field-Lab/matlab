figure
nds='87654321'
clear returns returnsWF
for i=1:8
    subplot(4,2,i)
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:249
        [~,m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end

    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:249
        [~,m]=max(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)<a(j),1);
        if ~isempty(b)
            returnsWF(j,i)=b+m;
        else
            returnsWF(j,i)=2000;
        end
    end
%     hist([returns(:,i)-returnsWF(:,i); -2000; 2000],80)
    
    hist([returns(:,i); 0; 2000],40)
%     axis([0 2000 0 30])
%     plot(returns(:,i),returnsWF(:,i),'o')
%     xlabel('Black latency')
%     ylabel('White latency')
    title(['ND',nds(i),' black - white'])
end


figure
for i=1:8
    subplot(2,4,i)
    [a,b]=hist([returns(:,i);-250;2150],24);
    plot(b(2:end-1),a(2:end-1)/2.49,'linewidth',2)
    hold on
    [a,b]=hist([returnsWF(:,i);-250;2150],24);
    plot(b(2:end-1),a(2:end-1)/2.49,'r','linewidth',2)
    axis([0 2010 0 70])
    title(['ND',nds(i)])
    legend('pos','neg','location','northwest')
    ylabel('% of all ON cells')
    xlabel('latency')
end


figure
for i=1:8
    subplot(2,4,i)
    m=find(returns(:,i)<1000);
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp(:,m),2));
    hold on
    m1=find(returns(:,i)>1999);
    plot(mean(tmp(:,m1),2),'r');
    m2=find(returns(:,i)<500&returns(:,i)>300);
    plot(mean(tmp(:,m2),2),'m');
%     axis([0 2010 0 70])
    title(['ND',nds(i),' tr ',int2str(length(m)),'  sus ',int2str(length(m1))])
%     legend('neg','pos','location','northwest')
%     ylabel('% of all ON cells')
%     xlabel('latency')
end


%Fig R4
figure
for i=1:8
    subplot(2,4,i)
    m=find(returns(:,i)<300);
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp(:,m),2));
    hold on
    m1=find(returns(:,i)>1999);
    plot(mean(tmp(:,m1),2),'r');
    axis([0 4500 -15 65])
    title(['ND',nds(i),' trans ',int2str(length(m)),'  sust ',int2str(length(m1))])
    legend(['<300ms'],['2000+ms'])
    ylabel('Hz')
    xlabel('time')
    line([1 500],[50,50]-10,'color','k')
    line([500 2500],[55,55]-10,'color','k')
    line([2500 4500],[50,50]-10,'color','k')
    line([500 500],[55,50]-10,'color','k')
    line([2500 2500],[55,50]-10,'color','k')
end



% Fig R3

figure
nds='87654321'
for i=1:8
    subplot(4,2,i)
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [k(i,j),m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end
    
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [kWF(i,j),m]=max(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)<a(j),1);
        if ~isempty(b)
            returnsWF(j,i)=b+m;
        else
            returnsWF(j,i)=2000;
        end
    end
    medWF(i)=mean(returnsWF(returnsWF(:,i)<1000,i));
    medBF(i)=mean(returns(returns(:,i)<1000,i));
    
    hist([returns(:,i); 0; 2000],40)
    axis([0 2000 0 100])
    %     plot(returns(:,i),returnsWF(:,i),'o')
    %     xlabel('Black latency')
    %     ylabel('White latency')
    title(['ND',nds(i),' black'])
end


figure
errorbar(mean(k'),std(k')/sqrt(200))
hold on
errorbar(mean(kWF'),std(kWF')/sqrt(200),'r')
for i=1:8
    [h,p(i)]=ttest(k(i,:),kWF(i,:));
    if h
        plot(i, mean(kWF(i,:))+10,'*k')
    end
end
xlabel('ND')
ylabel('HZ')
title('Mean +/-st error of peak amplitude of ON response across population')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
legend('Positive step','Negative step')


figure
plot(returns(:,3:5)')







%Fig R5
figure
col='bgr'
clear p
subtrspont=1;
for j=1:3
    for i=3:5
        subplot(1,3,j)
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
        %     sus=returns(:,i)>1999&kf&(~ff)&(~inh);
        m0=find(trans);
        m=find(inh);
        m1=find(sus);
        if ~subtrspont
            tmp=wf{i}(:,onOff>0);
        end
        switch j
            case 1
                tit='transient cells';
                toPlot=mean(tmp(:,m0),2);   
                p(i-2)=length(m0)/2.49;
            case 2
                tit='inhibitory cells';
                toPlot=mean(tmp(:,m),2);   
                p(i-2)=length(m)/2.49;
            case 3
                tit='sustained cells';
                toPlot=mean(tmp(:,m1),2);   
                p(i-2)=length(m1)/2.49;
        end
        
        plot(toPlot,col(i-2),'Linewidth',2);
        hold on
        axis([0 5000 -20 60])
        title(tit)
        
    end
    legend('ND6','ND5','ND4')
    for i=1:3
        bar(4500+i*130,p(i),120,'facecolor',col(i))
    end
end





% amount of sustained cells: ON, OFF
figure
for i=1:8
    subplot(2,4,i)

    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    kf=mean(tmp(2000:2450,:));
    kf=find(kf'>5);
    hold on
    plot(mean(tmp(:,kf),2),'r','linewidth',2)    
    kf1=kf;

    axis([0 4850 -20 65])
    title(['ND',nds(i),' sust ',int2str(round(length(kf)/2.49)),'%'])
    
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    kf=mean(tmp(2000:2450,:));
    kf=find(kf'>5);
    hold on
    plot(mean(tmp(:,kf),2),'color','b','linewidth',2) 
    legend('ON','OFF')

    axis([0 4850 -20 65])
    title(['ND',nds(i),' ON ',int2str(round(length(kf1)/2)),'% OFF ',int2str(round(length(kf)/2)),'%'])
    
    ylabel('Hz')
    xlabel('time')

    line([500 2500],[61,61],'color','k')

end


% Fig.R7
cellTypes=zeros(3,249);
for i=3:5
    
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
    
    cellTypes(i-2,m0)=1;
    cellTypes(i-2,m)=2;
    cellTypes(i-2,m1)=3;

end
for i=1:3
    stableCells(i)=sum(sum(cellTypes==i)==3)/2.49;
    almoststableCells(i)=sum(sum(cellTypes==i)==2)/2.49;
end
figure
bar([stableCells; almoststableCells]')
legend('same type at 3 NDs','same type at 2 NDs','location','northwest')
set(gca,'xticklabel',{'transient','inhibitory','sustained'})
ylabel('% of all ON cells')
title('% of ON cells maintaining the same type of response at ND6-4')

% sustained+inhibitory cells statistics
cellTypes=zeros(3,249);
for i=3:5
    
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
    
    cellTypes(i-2,m0)=1;
    cellTypes(i-2,m)=2;
    cellTypes(i-2,m1)=2;

end
clear stableCells almoststableCells
for i=1:2
    stableCells(i)=sum(sum(cellTypes==i)==3)/2.49;
    almoststableCells(i)=sum(sum(cellTypes==i)==2)/2.49;
end

%% OFF response

% Fig R8
figure
bord=1;
for i=1:8
    subplot(2,4,i)
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    
    ff=mean(tmp(2550:3000,:))';  
    tmp=wf{i}(:,onOff>0);
    hold on
    plot(mean(tmp(:,ff<-bord),2),'color',[0.3 0.8 0.1],'linewidth',2)
    plot(mean(tmp(:,ff>bord),2),'linewidth',2)
    plot(mean(tmp(:,ff>=-bord&ff<=bord),2),'r','linewidth',2)
    bar(4550,sum(ff<-bord)/2.49,90,'facecolor',[0.3 0.8 0.1])
    bar(4650,sum(ff>bord)/2.49,90,'facecolor','b')
    bar(4750,sum(ff>=-bord&ff<=bord)/2.49,90,'facecolor','r')

    
    axis([0 4850 -20 80])
    title(['ND',nds(i)])
    legend('<-1','>-1&<1','>1','location','northwest')
    ylabel('Hz')
    xlabel('time')
    line([1 500],[50,50]-5,'color','k')
    line([500 2500],[55,55]-5,'color','k')
    line([2500 4500],[50,50]-5,'color','k')
    line([500 500],[55,50]-5,'color','k')
    line([2500 2500],[55,50]-5,'color','k')

end

%% Asymmetry
%Fig R9
figure
bord=1;
% for end of the step, take 2000:2450 for ff
for i=1:8
    subplot(2,4,i)
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    tmp1=bf{i}(:,onOff>0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);

    
    ff=mean(tmp1(550:2000,:))';  
    tmp1=bf{i}(:,onOff>0);
    
    hold on
    plot(mean(tmp1(:,ff<-bord),2),'color',[0.3 0.8 0.1],'linewidth',2)
    plot(mean(tmp1(:,ff>bord),2),'linewidth',2)
    plot(mean(tmp1(:,ff>=-bord&ff<=bord),2),'r','linewidth',2)
    bar(4550,sum(ff<-bord)/2.49,90,'facecolor',[0.3 0.8 0.1])
    bar(4650,sum(ff>bord)/2.49,90,'facecolor','b')
    bar(4750,sum(ff>=-bord&ff<=bord)/2.49,90,'facecolor','r')

    
    axis([0 4850 -20 90])
    title(['ND',nds(i)])
    legend('<-1','>-1&<1','>1','location','northwest')
    ylabel('Hz')
    xlabel('time')
    line([1 500],[50,50]-5,'color','k')
    line([500 2500],[45,45]-5,'color','k')
    line([2500 4500],[50,50]-5,'color','k')
    line([500 500],[45,50]-5,'color','k')
    line([2500 2500],[45,50]-5,'color','k')

end

% Fig R10 - asymmetry (R6 and black flash)
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
%     sus=returns(:,i)>1999&kf&(~ff)&(~inh);
    m0=find(trans);
    m=find(inh);
    m1=find(sus);
    
    hold on
    bar(100,length(m0)/2.49,80,'facecolor',[0.3 0.8 0.1])
    bar(200,length(m)/2.49,80,'facecolor','b')
    bar(300,length(m1)/2.49,80,'facecolor','r')
    
    

    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    
    ff=mean(tmp((800:2450)+2000,:))';
    ff=ff<1;
    kf=mean(tmp((500:700)+2000,:))';
    fi=mean(tmp((680:820)+2000,:))';
    ti=mean(tmp((880:1000)+2000,:))';
    si=mean(tmp((2000:2450)+2000,:))';
    si=si>1;
    ti=ti>fi;
    fi=fi<(kf-3);
    kf=kf>3;
    trans=kf&ff;
    inh=fi&kf&(~ff)&ti;
    sus=si&kf&(~ff)&(~inh);
%     sus=returns(:,i)>1999&kf&(~ff)&(~inh);
    m0=find(trans);
    m=find(inh);
    m1=find(sus);   

    bar(700,length(m0)/2.49,80,'facecolor',[0.3 0.8 0.1])
    bar(800,length(m)/2.49,80,'facecolor','b')
    bar(900,length(m1)/2.49,80,'facecolor','r')
    ylabel('%, of all ON cells')
    set(gca,'xtick',200:600:800,'xticklabel',{'White Step','Black Step'})
    
    axis([0 1100 0 70])
    title(['ND',nds(i)])
    legend('trans','inh','sust')
end




figure
for i=1:8
    subplot(2,4,i)
    
    tmp=wf{i}(:,onOff>0);

    pl=mean(tmp(50:450,:))';    
    ff=max(tmp(500:900,:))';
    kf=min(tmp(2500:2900,:))';
    hist((kf-pl)./((kf+pl)))
    hold on
%     plot(sort(kf),'r')
%     
%     tmp=bf{i}(:,onOff<0);
%     tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
%     
%     ff=sum(tmp(500:900,:))';
%     kf=sum(tmp(2500:2900,:))';
%      plot(kf./ff,'.r')

    title(['ND',nds(i)])

end





%% trying to find groups of ON cells
% constant cells
cums=[];
cc=1;
types=cell(249,1);
for j=1:249
    tmp=cellTypes(:,j);
    cums=[cums j];
    for k=j+1:249
        if isempty(find(cums==k))&&sum(cellTypes(:,k)==tmp)==3
            types{cc}=[types{cc} k];
            cums=[cums,k];
        end
    end
    cc=cc+1;
end
cc=[];
for j=1:249
    if length(types{j})>10
        cc=[cc j];
        length(types{j})
    end
end

figure
lm=1;
for j=cc
    cellTypes(:,j)'
    for i=3:5
        subplot(5,3,lm)
        tmp=wf{i}(:,onOff>0);
        tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
        tmp=tmp(:,types{j});
        
        plot(mean(tmp,2),col(i-2))
        hold on
        
        lm=lm+1;
    end
end
figure
imagesc(cellTypes)

i=2
a=cellTypes(2,:)==i;
b=cellTypes(3,:)==i;
sum(a&b)/2.49
[a,b]=sort(cellTypes(1,:));
figure
imagesc(cellTypes(:,b))

figure
bord=1;
clear binind
for i=3:5
    tmp=wf{i}(:,onOff>0);
%     tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    if i==3
        r1=680:750;
        r2=900:1200;
    elseif i==4
        r1=670:900;
        r2=1000:1700;
    else
        r1=670:870;
        r2=900:1400;        
    end
    tt=mean(tmp(r1,:));
    mm=mean(tmp(r2,:));
    a=mm-tt;
    plot(sort(a),col(i-2))
    hold on
    line([0 249],[0,0],'color','k')
    line([0 249],[2,2],'color','k')

end



figure
bord=1;
clear binind
for i=3:5
    subplot(1,3,i-2)
    tmp=wf{i}(:,onOff>0);
%     tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    if i==3
        r1=680:750;
        r2=900:1200;
    elseif i==4
        r1=670:900;
        r2=1000:1700;
    else
        r1=670:870;
        r2=900:1400;        
    end
    tt=mean(tmp(r1,:));
    mm=mean(tmp(r2,:));
    a=mm-tt;
    plot(mean(tmp(:,a<-bord),2))
    p(i,1)=sum(a<-bord);
    p(i,2)=sum(a>bord);
    hold on
    plot(mean(tmp(:,a>bord),2),'r')
    binind(i-2,a>bord)=1;
    binind(i-2,a<-bord)=0;
    binind(i-2,a>=-bord&a<=bord)=8;

end

binind=binind(2:end,:)';
rt=sum(binind')<8;
binind=int2str(binind(rt,:));
binind=binind(:,1:3:end);
a=bin2dec(binind);

figure
clear means
cc=1;tt=1;
for code=0:3
    tk=a==code;
    sum(tk);
    if sum(tk)>10
        cc=tt;
        for i=3:5
            if i==3
                code
                exx=exp_codes(onOff>0);
                exx=exx(tk)
            end
            subplot(3,4,cc)
            kwf=wf{i}(:,onOff>0);
            kwf=kwf(:,rt);
            kwf=kwf-repmat(mean(kwf(50:450,:)),4500,1);            
            plot(kwf(:,tk))
            axis tight
            hold on
            plot(nanmean(kwf(:,tk),2),'color','k','linewidth',2)
            title([dec2bin(code,2),'   ',int2str(sum(tk)),' ',int2str(length(unique(exx)))])            
            cc=cc+4;
        end
        tt=tt+1;
    end
end




sum(binind(1,:)==1)



figure
i=5;
tmp=wf{i}(:,onOff>0);
r1=670:870;
r2=900:1400;
tt=mean(tmp(r1,:));
mm=mean(tmp(r2,:));
a=mm-tt;    
b=a>3;
plot(mean(tmp(:,b),2),'r')
hold on
tmp=wf{3}(:,onOff>0);
plot(mean(tmp(:,b),2),'b')
tmp=wf{4}(:,onOff>0);
plot(mean(tmp(:,b),2),'g')

figure
i=3;
tmp=wf{i}(:,onOff>0);
r1=670:870;
r2=900:1400;
tt=mean(tmp(r1,:));
mm=mean(tmp(r2,:));
a=mm-tt;    
b=a<-3;
plot(mean(tmp(:,b),2),'b')
hold on
tmp=wf{5}(:,onOff>0);
plot(mean(tmp(:,b),2),'r')
tmp=wf{4}(:,onOff>0);
plot(mean(tmp(:,b),2),'g')





figure
for i=3:5
    subplot(1,3,i-2)
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    tmp=sum(tmp(501:1000,:))./sum(tmp(501:2000,:));
    hist(tmp(tmp>0&tmp<1.3));

end


