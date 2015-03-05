tic
wt=[format_data('20120329',4,1); format_data('20120627',6,2);format_data('20120714',6,3);...
    format_data('20120902_1',4,4);format_data('20120902_2',4,5);format_data('20121023',12,6)];

cone=[format_data('20120927',12,1); format_data('20120827_1',4,2); format_data('20120828',4,3);format_data('20120920',12,4);...
    format_data('20120920_1',12,5);format_data('20120921',12,6);format_data('20120925',12,7)];

rod=[format_data('20120928',12,1);format_data('20121002',12,2)];

cpfl=format_data('20121018',12,1);

mela=[format_data('20120829_1',4,1);format_data('20120829_2',4,2);format_data('20120829_3',4,3);format_data('20120830',6,4)];

save('S:\Paper factory\Alex\2_lightAdaptation\data\collect_param.mat','wt','cone','rod','cpfl','mela')
toc

ms=20;

figure
set(gcf,'position',[1 31 1680 946])
for datType=1:3
    if datType==1
        pulled=wt;
        col='b';
    elseif datType==2
        pulled=rod;
        col='r';
    elseif datType==3
        pulled=cone;
        col='g';
    elseif datType==4
        pulled=cpfl;
        col='m';
    elseif datType==5
        pulled=mela;
        col='c';
    end
    subplot(3,1,1)
    hold on
    set(gca,'XTick',6:12:96,'xticklabel',{'8','7','6','5','4','3','2','1'})
    
    a=pulled(:,:,2);
    clear k
    for i=1:size(a,2)
        k(i)=mean(a(a(:,i)>0,i));
    end
    
    plot(k,'.','MarkerSize',ms,'color',col)
    
    for i=12.5:12:96
        line([i,i],[0, 200],'color','k')
    end
    axis([0,96,50,200])
    title('Latency')
    
    subplot(3,1,2)
    hold on
    set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
    a=pulled(:,:,1);
    clear k n
    for i=1:size(a,2)
        k(i)=mean(a(a(:,i)>0,i));
        n(i)=mean(a(a(:,i)<0,i));
    end
    
    plot(k,'.','MarkerSize',ms,'color',col)
    hold on
    plot(n,'.','MarkerSize',ms,'color',col)
    
    for i=12.5:12:146
        line([i,i],[-400, 400],'color','k')
    end
    line([0,96],[0,0],'color','k')
    axis([0,96,-400,400])
    title('Amplitude')
    
    
    subplot(3,1,3)
    hold on
    set(gca,'XTick',6:12:96,'xticklabel',{'8','7','6','5','4','3','2','1'})
    
    a=pulled(:,:,3);
    clear k
    for i=1:size(a,2)
        k(i)=mean(a(a(:,i)>0,i));
    end
    
    plot(k,'.','MarkerSize',ms,'color',col)
    
    for i=12.5:12:96
        line([i,i],[0, 60],'color','k')
    end
    axis([0,96,0,60])
    title('Width')
    
end


% select parameter
par=9;
sel=6;
ms=20;

figure
set(gcf,'position',[1 31 1680 946])
for datType=1:5
    if datType==1
        unitSelect=sum(wt(:,:,par)==sel,2)>0;
        pulled=wt(unitSelect,:,:);
        col='b';
    elseif datType==2
        unitSelect=sum(rod(:,:,par)==sel,2)>0;
        pulled=rod(unitSelect,:,:);
        col='r';
    elseif datType==3
        unitSelect=sum(cone(:,:,par)==sel,2)>0;
        pulled=cone(unitSelect,:,:);
        col='g';
    elseif datType==4
        unitSelect=sum(cpfl(:,:,par)==sel,2)>0;
        pulled=cpfl(unitSelect,:,:);
        col='m';
    elseif datType==5
        unitSelect=sum(mela(:,:,par)==sel,2)>0;
        pulled=mela(unitSelect,:,:);
        col='c';
    end
    subplot(3,1,1)
    hold on
    set(gca,'XTick',6:12:96,'xticklabel',{'8','7','6','5','4','3','2','1'})
    
    a=pulled(:,:,2);
    clear k
    for i=1:size(a,2)
        k(i)=mean(a(a(:,i)>0,i));
    end
    
    plot(k,'.','MarkerSize',ms,'color',col)
    
    for i=12.5:12:96
        line([i,i],[0, 200],'color','k')
    end
    axis([0,96,50,200])
    title('Latency')
    
    subplot(3,1,2)
    hold on
    set(gca,'XTick',6:12:146,'xticklabel',{'8','7','6','5','4','3','2','1','2','3','4','5'})
    a=pulled(:,:,1);
    clear k n
    for i=1:size(a,2)
        k(i)=mean(a(a(:,i)>0,i));
        n(i)=mean(a(a(:,i)<0,i));
    end
    
    plot(k,'.','MarkerSize',ms,'color',col)
    hold on
    plot(n,'.','MarkerSize',ms,'color',col)
    
    for i=12.5:12:146
        line([i,i],[-400, 400],'color','k')
    end
    line([0,96],[0,0],'color','k')
    axis([0,96,-400,400])
    title('Amplitude')
    
    
    subplot(3,1,3)
    hold on
    set(gca,'XTick',6:12:96,'xticklabel',{'8','7','6','5','4','3','2','1'})
    
    a=pulled(:,:,3);
    clear k
    for i=1:size(a,2)
        k(i)=mean(a(a(:,i)>0,i));
    end
    
    plot(k,'.','MarkerSize',ms,'color',col)
    
    for i=12.5:12:96
        line([i,i],[0, 60],'color','k')
    end
    axis([0,96,0,60])
    title('Width')
    
end



AxisWTConeRod=[1.5,8.5,0,45; 1.5,8.5,-15,25; 1.5,8.5,-10,120];
% differences
par=9;
sel=12;
ms=20;
tr=[2,12;1,11;1,12];

figure
set(gcf,'position',[1 31 1680 946])

for sbpl=1:3
    subplot(1,3,sbpl)
    if sbpl==1
        tit='WT';
        pulled=wt;
    elseif sbpl==2
        tit='CONE';
        pulled=cone;
    else
        tit='ROD';
        pulled=rod;
    end
    clear k p
    b=cell(3,8);
    cnt=1;
    for sel=[4 6 12]
        unitSelect=sum(pulled(:,:,par)==sel,2)>0;
        if sum(unitSelect)>0
            a=pulled(unitSelect,:,2);
            c=a(:,tr(cnt,2):12:end);
            d=a(:,tr(cnt,1):12:end);
            %     c=(a(:,11:12:end)+a(:,12:12:end))/2;
            %     d=(a(:,1:12:end)+a(:,2:12:end))/2;
            a=c-d;
            
            for i=1:8
                b{cnt,i}=a(~isnan(a(:,i)),i);
            end
        end
        cnt=cnt+1;
    end
    clear k p
    for i=1:8
        a=[b{1,i}; b{2,i}; b{3,i}];
        m(i)=numel(a);
        k(i,1)=mean(a);
        k(i,2)=std(a)/numel(a);
        [~,p(i)]=ttest(a);
    end
    bar(k(:,1));
    hold on
    errorbar(k(:,1),k(:,2),'xr','LineWidth',2)
    
    for i=2:8
        text(i-0.5,k(i,1)+sign(k(i,1))*(k(i,2)+4),num2str(p(i)))
        text(i-0.1,k(i,1)+sign(k(i,1))*(k(i,2)+2),num2str(m(i)))
    end
    set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
    axis(AxisWTConeRod(sbpl,:))
    title(tit,'fontsize',16,'fontweight','b')
end






AxisWTConeRod=[1.5,8.5,-90,25; 1.5,8.5,-15,25; 1.5,8.5,-10,120];
% differences
par=9;
sel=12;
ms=20;
tr=[2,12;1,12];

figure
set(gcf,'position',[1 31 1680 946])
sel=[4 12]
for sbpl=1:2
    subplot(1,2,sbpl)
    tit='CONE';
    pulled=cone;
    clear k p
    b=cell(1,8);
    unitSelect=sum(pulled(:,:,par)==sel(sbpl),2)>0;
    if sum(unitSelect)>0
        a=pulled(unitSelect,:,2);
        c=a(:,tr(sbpl,2):12:end);
        d=a(:,tr(sbpl,1):12:end);
        %     c=(a(:,11:12:end)+a(:,12:12:end))/2;
        %     d=(a(:,1:12:end)+a(:,2:12:end))/2;
        a=c-d;
        
        for i=1:8
            b{i}=a(~isnan(a(:,i)),i);
        end
        
        clear k p
        for i=1:8
            a=b{1,i};
            m(i)=numel(a);
            k(i,1)=mean(a);
            k(i,2)=std(a)/numel(a);
            [~,p(i)]=ttest(a);
        end
        bar(k(:,1));
        hold on
        errorbar(k(:,1),k(:,2),'xr','LineWidth',2)
        
        for i=2:8
            text(i-0.5,k(i,1)+sign(k(i,1))*(k(i,2)+4),num2str(p(i)))
            text(i-0.1,k(i,1)+sign(k(i,1))*(k(i,2)+2),num2str(m(i)))
        end
        set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
        axis(AxisWTConeRod(sbpl,:))
        title([tit, '  trials: ',int2str(sel(sbpl))],'fontsize',16,'fontweight','b')
    end
end






AxisWTConeRod=[1.5,8.5,-90,25; 1.5,8.5,-15,25; 1.5,8.5,-10,120];
% differences
par=9;
sel=12;
ms=20;
tr=[2,12;1,12];

figure
set(gcf,'position',[1 31 1680 946])
sel=[4 12]
for sbpl=1:2
    subplot(1,2,sbpl)
    tit='CONE';
    pulled=cone;
    clear k p
    b=cell(1,8);
    unitSelect=sum(pulled(:,:,par)==sel(sbpl),2)>0;
    if sum(unitSelect)>0
        a=pulled(unitSelect,:,2);
        c=a(:,tr(sbpl,2):12:end);
        d=a(:,tr(sbpl,1):12:end);
        %     c=(a(:,11:12:end)+a(:,12:12:end))/2;
        %     d=(a(:,1:12:end)+a(:,2:12:end))/2;
        a=c-d;
        
        for i=1:8
            b{i}=a(~isnan(a(:,i)),i);
        end
        
        clear k p
        for i=1:8
            a=b{1,i};
            m(i)=numel(a);
            k(i,1)=mean(a);
            k(i,2)=std(a)/numel(a);
            [~,p(i)]=ttest(a);
        end
        bar(k(:,1));
        hold on
        errorbar(k(:,1),k(:,2),'xr','LineWidth',2)
        
        for i=2:8
            text(i-0.5,k(i,1)+sign(k(i,1))*(k(i,2)+4),num2str(p(i)))
            text(i-0.1,k(i,1)+sign(k(i,1))*(k(i,2)+2),num2str(m(i)))
        end
        set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
        axis(AxisWTConeRod(sbpl,:))
        title([tit, '  trials: ',int2str(sel(sbpl))],'fontsize',16,'fontweight','b')
    end
end



figure;plot(b{6},'.')



AxisWTConeRod=[1.5,8.5,-100,25; 1.5,8.5,-100,25; 1.5,8.5,-10,120];
% differences
par=9;
sel=12;
ms=20;
tr=[2,12;1,12];

figure
set(gcf,'position',[1 31 1680 946])
sel=[4 12]
for sbpl=1:2
    subplot(1,2,sbpl)
    tit='CONE';
    pulled=cone;
    clear k p
    b=cell(1,8);
    unitSelect=sum(pulled(:,:,par)==sel(sbpl),2)>0;
    if sum(unitSelect)>0
        a=pulled(unitSelect,:,1);
        c=a(:,tr(sbpl,2):12:end);
        d=a(:,tr(sbpl,1):12:end);
        %     c=(a(:,11:12:end)+a(:,12:12:end))/2;
        %     d=(a(:,1:12:end)+a(:,2:12:end))/2;
        a=(abs(c)-abs(d))./abs(d);
        
        for i=1:8
            b{i}=a(~isnan(a(:,i)),i);
        end
        
        clear k p
        for i=1:8
            a=b{1,i};
            m(i)=numel(a);
            k(i,1)=mean(a);
            k(i,2)=std(a)/numel(a);
            [~,p(i)]=ttest(a);
        end
        bar(k(:,1));
        hold on
        errorbar(k(:,1),k(:,2),'xr','LineWidth',2)
        
        for i=2:8
            text(i-0.5,k(i,1)+sign(k(i,1))*(k(i,2)+4),num2str(p(i)))
            text(i-0.1,k(i,1)+sign(k(i,1))*(k(i,2)+2),num2str(m(i)))
        end
        set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
        axis(AxisWTConeRod(sbpl,:))
        title([tit, '  trials: ',int2str(sel(sbpl))],'fontsize',16,'fontweight','b')
    end
end







AxisWTConeRod=[1.5,8.5,0,45; 1.5,8.5,-15,25; 1.5,8.5,-10,120];
% differences
par=9;
sel=12;
ms=20;
tr=[2,12;1,11;1,12];

figure
set(gcf,'position',[1 31 1680 946])

for sbpl=1:3
    subplot(1,3,sbpl)
    if sbpl==1
        tit='WT';
        pulled=wt;
    elseif sbpl==2
        tit='CONE';
        pulled=cone;
    else
        tit='ROD';
        pulled=rod;
    end
    clear k p
    b=cell(3,8);
    cnt=1;
    for sel=[4 6 12]
        unitSelect=sum(pulled(:,:,par)==sel,2)>0;
        if sum(unitSelect)>0
            a=pulled(unitSelect,:,1);
            c=a(:,tr(cnt,2):12:end);
            d=a(:,tr(cnt,1):12:end);
            %     c=(a(:,11:12:end)+a(:,12:12:end))/2;
            %     d=(a(:,1:12:end)+a(:,2:12:end))/2;
            a=(abs(c)-abs(d))./abs(d);
            
            for i=1:8
                b{cnt,i}=a(~isnan(a(:,i)),i);
            end
        end
        cnt=cnt+1;
    end
    clear k p
    for i=1:8
        a=[b{1,i}; b{2,i}; b{3,i}];
        m(i)=numel(a);
        k(i,1)=mean(a);
        k(i,2)=std(a)/numel(a);
        [~,p(i)]=ttest(a);
    end
    bar(k(:,1));
    hold on
    errorbar(k(:,1),k(:,2),'xr','LineWidth',2)
    
    for i=2:8
        text(i-0.5,k(i,1)+sign(k(i,1))*(k(i,2)+4),num2str(p(i)))
        text(i-0.1,k(i,1)+sign(k(i,1))*(k(i,2)+2),num2str(m(i)))
    end
    set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
    axis(AxisWTConeRod(sbpl,:))
    title(tit,'fontsize',16,'fontweight','b')
end
