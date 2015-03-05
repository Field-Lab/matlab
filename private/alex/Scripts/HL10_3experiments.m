cd('/mnt/muench_data/user/alexandra/scripts')

clear
codeWord='HL10'
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])

%% Calculate parameters


% choose cells with stable firing rate (by std), for each ND separately
% meaningful value: 25
cellByND=zeros(16,103);
figure
for k=1:16
    subplot(4,4,k)
    hold on
    for cellID=1:103
        t=reshape(SpikeCount(k,cellID,:),dim,1);
        plot(cellID,std(t(1:2:end)),'*','color',col(cellID));
        if std(t(1:2:end))<25
            cellByND(k,cellID)=cellID;
        end
    end
    line([0 104],[25,25],'color','k')
    axis([0 104 0 50])
    title(['trial',int2str(k)])
end


correlFiringRate=zeros(16,103); % firing rate vs contrast
correlIndex=zeros(16,103); % zero crossing vs contrast
correlNonlin=zeros(16,103); % gain vs contrast
spikeCount_cum=zeros(16,103); % summed firing rate (NB! exceeds real value due to overlap of individual intervals)
correlNonlinAdj=zeros(16,103); % gain divided by filter amplitude (peak value) vs contrast
for k=1:16
    for cellID=1:103
        t=reshape(SpikeCount(k,cellID,:),dim,1);
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

%% PLOTS


figure
nums=zeros(16,2);
nds='8877665544332211';
for k=1:16
    subplot(4,4,k)
    hold on
    a=correlFiringRate(k,onOff>0&cellByND(k,:)>0);
%     a=correlContrast(k,onOff>0);
    plot(a,'r*')
    errorbar(58,mean(a(a>0.8)),std(a(a>0.8)),'.r','MarkerSize',20)
    errorbar(58,mean(a(a<-0.8)),std(a(a<-0.8)),'.r','MarkerSize',20)
    nums(k,1)=sum(a<-0.8);
    nums(k,2)=numel(a);

    b=correlFiringRate(k,onOff<0&cellByND(k,:)>0);
%     b=correlContrast(k,onOff<0);
    plot(b,'b*')    
    errorbar(62,mean(b(b>0.8)),std(b(b>0.8)),'.','MarkerSize',20)
    errorbar(62,mean(b(b<-0.8)),std(b(b<-0.8)),'.','MarkerSize',20)
    
    axis([0 65 -1.2 1.2])
    line([0 65],[0,0],'color','k')
    line([0 65],[0.8,0.8],'color','k')
    line([0 65],[-0.8,-0.8],'color','k')
    title(['ND', nds(k),', red ON, blue OFF'])
end




for k=1:10000
    j(k)=corr(sigALL',rand(6,1));
end

figure
plot(j,'*')
hist(j,50)
