cd('/mnt/muench_data/user/alexandra/scripts')

clear
codeWord='sine'
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])


%% Calculate parameters



for k=1:8
    t=reshape(ind(:,k,:),dim,size(LinearFilter,4) );
    for j=1:dim
        meanInd(j,k)=mean(t(j,t(j,:)>70&t(j,:)<250));
    end
end
figure
plot(sigALL,meanInd)
legend(['7654321']')



% correlation of parameters with contrast changes
correlFiringRate=zeros(8,103); % firing rate vs contrast
correlIndex=zeros(8,103); % zero crossing vs contrast
correlNonlin=zeros(8,103); % gain vs contrast
spikeCount_cum=zeros(8,103); % summed firing rate (NB! exceeds real value due to overlap of individual intervals)
correlNonlinAdj=zeros(8,103); % gain divided by filter amplitude (peak value) vs contrast
for k=1:8
    for cellID=1:103
        t=reshape(SpikeCount(k,cellID,:),272,1);
        correlFiringRate(k,cellID)=corr(sigALL',t);        
        spikeCount_cum(k,cellID)=sum(t);
        t=ind(:,k,cellID);
        correlIndex(k,cellID)=corr(sigALL',t);
        t=nonLinear(2,:,k,cellID)';
        correlNonlin(k,cellID)=corr(sigALL',t);
        t=nonLinear_adjusted(:,k,cellID);
        correlNonlinAdj(k,cellID)=corr(sigALL',t);
    end
end


% choose cells with stable firing rate (by std), for each ND separately
% meaningful value: 15
cellByND=zeros(8,103);
periodsTake=58;
figure
for k=1:8
    subplot(2,4,k)
    hold on
    for cellID=1:103
        t=reshape(SpikeCount(k,cellID,:),272,1);
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
        plot(cellID,std(t(20:40)),'*','color',col(cellID));
        if std(t(20:40))<15
            cellByND(k,cellID)=cellID;
        end
    end
    line([0 104],[15,15],'color','k')
    title(['ND',int2str(9-k)])
end

% remove cells with unstable filters (high std) at high contrasts at ND6
% threshold value: 7

cellByND_ND6=cellByND;
ad=zeros(1,103);
figure

hold on
for cellID=1:103
    t=ind(:,3,cellID);
    t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
    t=mean(t');
    t=t-mean(t([1:3, 56:58]));
    hc=mean(t([1:3, 56:58]));
    lc=mean(t(23:33));
    
    if mean(t(23:33))>0
        ad(cellID)=1;
    end
    if onOff(cellID)>0
        cc='r';
    else
        cc='b';
    end
    
    if std(t([1:10, 48:58]))>7
        cellByND_ND6(:,cellID)=0;        
        %         plot(t)
    end
    plot(std(t([1:10, 48:58])),lc-hc,'*','color',cc)
    
    if onOff(cellID)>0
        col(cellID)='r';
    else
        col(cellID)='b';
    end
end
line([0 7],[0,14],'color','k')
line([0 7],[0,-14],'color','k')
line([0 20],[0,0],'color','k')
line([7 7],[-14,14],'color','k')

%% PLOTS


figure
for k=1:8
    subplot(2,4,k)
    hold on
    line([0 650],[0,0],'color','k')
    for cellID=find(cellByND_ND6(k,:))
        t=reshape(peak(41-13:end-13,k,cellID),periodsTake,floor(272/periodsTake));
        t(1:29,5:8)=t(end:-1:30,:);
        t(30:end,:)=[];  
        t=mean(t');
        p=reshape(SpikeCount(k,cellID,:),272,1);
        p=reshape(p(41-13:end-13),periodsTake,floor(272/periodsTake));
        p(1:29,5:8)=p(end:-1:30,:);
        p(30:end,:)=[];        
        p=mean(p');
        plot(p,t,col(cellID),'marker','.');
        if correlFiringRate(k,cellID)<-0.2;
            plot(p,t,'k','marker','.');
        end
    end
    axis([0 650 -100 100])
    title(['ND',int2str(9-k)])
    xlabel('spike count')
    ylabel('filter peak')
end


figure
nums=zeros(8,2);
for k=1:8
    subplot(2,4,k)
    hold on
    a=correlFiringRate(k,onOff>0&cellByND_ND6(k,:)>0);
%     a=correlContrast(k,onOff>0);
    plot(a,'r*')
    errorbar(58,mean(a(a>0.2)),std(a(a>0.2)),'.r','MarkerSize',20)
    errorbar(58,mean(a(a<-0.2)),std(a(a<-0.2)),'.r','MarkerSize',20)
    nums(k,1)=sum(a<-0.2);
    nums(k,2)=numel(a);

    b=correlFiringRate(k,onOff<0&cellByND_ND6(k,:)>0);
%     b=correlContrast(k,onOff<0);
    plot(b,'b*')    
    errorbar(62,mean(b(b>0.2)),std(b(b>0.2)),'.','MarkerSize',20)
    errorbar(62,mean(b(b<-0.2)),std(b(b<-0.2)),'.','MarkerSize',20)
    
    axis([0 65 -1.2 1.2])
    line([0 65],[0,0],'color','k')
    line([0 65],[0.2,0.2],'color','k')
    line([0 65],[-0.2,-0.2],'color','k')
    title(['ND', int2str(9-k),', red ON, blue OFF'])
end


figure
r=1;
for k=2:7
    for cellID=find(cellByND_ND6(k,:))
        t=reshape(SpikeCount(k,cellID,:),272,1);
%         if onOff(cellID)<0
            if corr(sigALL',t)>0.3
                subplot(3,4,r)
                plot(t,col(cellID))
                if correlIndex(k,cellID)<-0.5&&onOff(cellID)>0
                    plot(t,'k')
                end
                hold on
                title(['ND',int2str(9-k)])
            elseif corr(sigALL',t)<-0.3
                subplot(3,4,r+1)
                plot(t,col(cellID))
                if correlIndex(k,cellID)<-0.5&&onOff(cellID)>0
                    plot(t,'k')
                end
                hold on
                title(['ND',int2str(9-k)])
            end
%         end
    end
    r=r+2;
end


% adaptability index
figure
for k=1:8
    subplot(2,4,k)
    hold on
    a=correlFiringRate(k,onOff>0&cellByND(k,:)>0);
    b=correlIndex(k,onOff>0&cellByND(k,:)>0);
    plot(a,b,'ro')


    a=correlFiringRate(k,onOff<0&cellByND(k,:)>0);
    b=correlIndex(k,onOff<0&cellByND(k,:)>0);
    plot(a,b,'bo')
    
    axis([-1 1 -1 1])
    line([-1 1],[0,0],'color','k')
    line([0 0],[-1 1],'color','k')
    xlabel('firing rate')
    ylabel('latency')
%     line([0 65],[-0.2,-0.2],'color','k')
    title(['ND', int2str(9-k),', red ON, blue OFF'])
end


%% Linear part plot
periodsTake=58;
cell2take=[];
ad=zeros(1,103);
figure
for k=1:8
    subplot(2,4,k)
    hold on
    for cellID=1:103
        t=ind(:,k,cellID);
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
        t=t-mean(t([1:3, 56:58]));
        hc=mean(t([1:3, 56:58]));
        lc=mean(t(23:33));
        
        if mean(t(23:33))>0
            ad(cellID)=1;
        end
        if onOff(cellID)>0
            cc='r';
        else
            cc='b';
        end
        
        if std(t([1:10, 48:58]))<7
            cell2take=[cell2take cellID];
            %         plot(t)
        end
        plot(std(t([1:10, 48:58])),lc-hc,'*','color',cc)
        
        if onOff(cellID)>0
            col(cellID)='r';
        else
            col(cellID)='b';
        end
    end
    line([0 7],[0,14],'color','k')
    line([0 7],[0,-14],'color','k')
    line([0 20],[0,0],'color','k')
    line([7 7],[-14,14],'color','k')
end







periodsTake=58;
cell2take=[];
ad=zeros(1,103);
figure
hold on
kk=-0.25:1/206:0.25;
for cellID=1:103
    for k=1:8
        t=ind(:,k,cellID);
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
        t=t-mean(t([1:3, 56:58]));
        hc=mean(t([1:3, 56:58]));
        lc=mean(t(23:33));
        
        if mean(t(23:33))>0
            ad(cellID)=1;
        end
        
        if std(t([1:10, 48:58]))<7
            cell2take=[cell2take cellID];
            %         plot(t)
        end
        zt(k)=lc-hc;        
    end
    if onOff(cellID)>0
        cc='r';
    else
        cc='b';
    end
    plot([1:8]+kk(cellID),zt,'color',cc)
    plot([1:8]+kk(cellID),zt,'o','color',cc)
end
set(gca,'xtick',1:8,'xticklabel',[8:-1:1])
line([0,9],[0,0],'color','k')






figure
acc_t=zeros(8,58,2);
cnt_t=zeros(8,2);
for cellID=cell2take
    for k=1:8
        t=ind(:,k,cellID);
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
%         t=t-mean(t(1:10));
        subplot(2,8,k+ad(cellID)*8)
%         line([0 59],[10,10],'color','k')
%         line([28.5,28.5],[-50 150],'color','k')
        plot(t,col(cellID));
        hold on
%         axis([0 59 -50 150])
%         if t(28)<10
%             acc_t(k,:,ad(cellID)+1)=t+acc_t(k,:,ad(cellID)+1);
%             cnt_t(k,ad(cellID)+1)=cnt_t(k,ad(cellID)+1)+1;
%         end
    end
end
col(wrongCell)='k'

for i=2:7
    subplot(2,6,i-1)
    plot(acc_t(i,:,1)./cnt_t(i,1),'k','linewidth',5)
    title(int2str(cnt_t(i,1)))
    
    subplot(2,6,i+5)
    plot(acc_t(i,:,2)./cnt_t(i,2),'k','linewidth',5)
    title(int2str(cnt_t(i,2)))
end


%% Linear part plot not separated
figure
for cellID=cell2take
    for k=1:7
        t=ind(:,k,cellID);
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
        t=t-mean(t(1:5));
        subplot(2,4,k)
        line([0 59],[0,0],'color','k')
        plot(t,col(cellID));
        hold on
    end
end
for k=1:7
    subplot(2,4,k)
    axis tight
    title(['ND',int2str(8-k)])
end


%% nonLinear part plot 


figure
acc_off=[];
acc_on=[];
for cellID=unique(cell2take)
    
    for k=1:8
        t=nonLinear(2,:,k,cellID);        
%       t=reshape(t(41:end),periodsTake,floor(272/periodsTake));
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
        t=ones(1,length(t))./t;
%       t=circshift(t,[0,-18]);
        subplot(2,4,k)
        hold on
        line([0 59],[0,0],'color','k')        
        plot(t,col(cellID));
        if onOff(cellID)<0
            acc_off=[acc_off; t];
        else
            acc_on=[acc_on; t];
        end
        axis([0 59 0 3000])
        
        title(['ND',int2str(9-k)])
    end
end


figure
for k=1:8
   
    subplot(2,4,k)
    
        plot(median(acc_off(k:8:end,:)),'linewidth',3,'color','b');
         hold on
          plot(median(acc_on(k:8:end,:)),'linewidth',3,'color','r');
          axis([0 59 0 0.4])
end



%% gain vs zero crossing

acc_zc=zeros(7,58,103);
acc_g=zeros(7,58,103);
for cellID=cell2take
    for k=1:7
        t=ind(:,k,cellID);
        t=reshape(t(41:end),periodsTake,floor(272/periodsTake)); 
        t=mean(t');
        t=circshift(t,[0,-18]);
        t=t-mean(t(1:10));
        acc_zc(k,:,cellID)=t;        
        
        t=nonLinear(2,:,k+1,cellID);
        t=reshape(t(41:end),periodsTake,floor(272/periodsTake));
        t=mean(t');
        t=circshift(t,[0,-18]);
        acc_g(k,:,cellID)=t;         
    end
end


figure

for cellID=cell2take
    for k=1:7
        subplot(2,4,k)
        hold on
        plot(acc_zc(k,:,cellID),acc_g(k,:,cellID),'.','color',col(cellID));
        title(['ND',int2str(8-k)])
    end
end



%% firing rate

figure
acc_off=[];
acc_on=[];
% cell2take_corr=cell2take;
% wrongCell=[];
z=zeros(8,2);
on=sum(onOff(cell2take_corr)==1);
off=sum(onOff(cell2take_corr)==-1);
for cellID=cell2take_corr
    
    for k=1:8
        t=reshape(SpikeCount(k,cellID,:),272,1);
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
%         t=circshift(t,[0,-18]);
        subplot(2,4,k)
        hold on
        line([0 59],[0,0],'color','k')
        
        tt=t-min(t);
        plot(tt,col(cellID));
        
        
        if (k==4&&max(tt)>150) || (k==3&&max(tt)>300) || (k==4&&tt(22)>100)
            cell2take_corr(cell2take_corr==cellID)=[];
        elseif (mean(tt(26:30))-mean(tt([1 2 57 58])))>10
                if onOff(cellID)>0
                    plot(tt,'k','LineWidth',2);
                    z(k,1)=z(k,1)+1;
                    if k==5
                        cellID
                    end
                else
                    plot(tt,'g','LineWidth',2);
                    z(k,2)=z(k,2)+1;
                end
%             wrongCell=[wrongCell cellID];
        end
        
        
        if onOff(cellID)<0
            acc_off=[acc_off; t];
        else
            acc_on=[acc_on; t];
        end
        axis([0 59 0 300])
        
        title(['ND',int2str(9-k)])
    end
end

wrongCell=unique(wrongCell);
col(wrongCell)='k'



%% firing rate per ND


        

figure

z=zeros(8,2);
firingRateChange=zeros(8,103);
for k=1:8
    subplot(2,4,k)
    hold on
    line([0 59],[0,0],'color','k')
    for cellID=cellByND(k,cellByND(k,:)>0)
        t=reshape(SpikeCount(k,cellID,:),272,1);
        t=reshape(t(41-13:end-13),periodsTake,floor(272/periodsTake));
        t=mean(t');
        tt=t-min(t);
        plot(tt,col(cellID));
        firingRateChange(k,cellID)=mean(tt([1:3 56:58]))-mean(tt(26:30));
        if (mean(tt(26:30))/mean(tt(1:3)))>1.5&&(mean(tt(26:30))/mean(tt(56:58)))>1.5
                if onOff(cellID)>0
                    plot(tt,'k','LineWidth',2);
                    z(k,1)=z(k,1)+1;
                else
                    plot(tt,'g','LineWidth',2);
                    z(k,2)=z(k,2)+1;
                end
        end
    end
    axis([0 59 0 300])    
    title(['ND',int2str(9-k)])
end

for k=1:8
    FiringRateOn(k)=z(k,1)/sum(onOff(cellByND(k,cellByND(k,:)>0))>0);
    FiringRateOff(k)=z(k,2)/sum(onOff(cellByND(k,cellByND(k,:)>0))<0);
end

% firing rate change
figure
for k=1:8
    subplot(2,4,k)
    hold on
    a=find(firingRateChange(k,cellByND(k,cellByND(k,:)>0)));
    a1=a(onOff(a)>0);
    a2=firingRateChange(k,a1);
    a2(isinf(a2))=nan;
    plot(a2,'r*')
    plot(50,nanmedian(a2(a2>1)),'r.','markersize',30)

    a1=a(onOff(a)<0);
    a2=firingRateChange(k,a1);
    a2(isinf(a2))=nan;
    plot(a2,'b*')
    plot(50,nanmedian(a2(a2>1)),'b.','markersize',30)
    axis([0 53 -50 200])
    line([0 53],[0,0],'color','k')
end





% plot it
figure

for k=1:8
    subplot(2,4,k)
    hold on
    line([0 59],[0,0],'color','k')
    for cellID=cellByND(k,cellByND(k,:)>0)
        t=reshape(nonLinear(1,41-13:end-13,k,cellID),periodsTake,floor(272/periodsTake));
        t=mean(t');
%         t=ones(1,58)./t;
%         plot(cellID,nanmean(t([1:3, 56:58]))/nanmean(t(26:30)),col(cellID),'marker','*');
        plot(t,col(cellID));
    end
    axis([0 59 0 50])    
    title(['ND',int2str(9-k)])
end






















for k=1:10000
    j(k)=corr(sigALL',randn(272,1));
end

figure
plot(j,'*')
hist(j,50)