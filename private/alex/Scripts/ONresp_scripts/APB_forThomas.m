addpath(genpath('/Users/alexth/Desktop/Scripts'))


clear
load('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat')


onCells=[];
offCells=[];
for i=1:size(datarun_nd,2)
    if datarun_nd(i).polarity==1
        onCells=[onCells i];
    elseif datarun_nd(i).polarity==-1
        offCells=[offCells i];
    end
end

%%
%%%%%%%%%% OFF cells analysis %%%%%%%%%%%%%

%% common parameters

reversed=0.3; % threshold for reversibility in wash 
stable=0.3; % threshold for consistency within condition  - stabilizing

%% Plot for diversity showing - black flash

conditions={'control','apb','wash'};

figure
contr='black';
cnt=1;
for i=1:3
    for nd=[6 5 4]        
        subplot(3,3,cnt)
        [data,~]=goodCellsSelection(datarun_nd, offCells, nd, contr, conditions{i}, 'threshRev', reversed,'threshStable', stable);
        
        plot(data)
        hold on
        plot(mean(data,2),'k','linewidth',3)
        title(['ND',int2str(nd),', ',contr,' flash, ',conditions{i},', n=' int2str(size(data,2))])
        line([500 500],[-1.5 1.5],'color','r')
        line([2500 2500],[-1.5 1.5],'color','r')
        axis([0 4500 -1.5 1.5])
        
        cnt=cnt+1;
    end
end



%%% weird cell

trustOff=intersect(offCells,find(ifReversable(:,1)>threshRev));
figure
for k=1:3
    subplot(3,1,k)
    acc=[];
    for i=trustOff(33)'
        data=datarun_nd(i).black.nd6(:,k);
        data=data-mean(data(1:450));
        data=data/max(data(500:1000));
        acc=[acc data];
    end
    plot(acc)
    hold on
    plot(mean(acc,2),'k','linewidth',4)
end



%% Plot for diversity showing - white flash


conditions={'control','apb','wash'};

figure
contr='white';
cnt=1;
for i=1:3
    for nd=[6 5 4]        
        subplot(3,3,cnt)
        [data,~]=goodCellsSelection(datarun_nd, offCells, nd, contr, conditions{i}, 'threshRev', reversed,'threshStable', stable);
 
        plot(data)
        hold on
        plot(mean(data,2),'k','linewidth',3)
        title(['ND',int2str(nd),', ',contr,' flash, ',conditions{i},', n=' int2str(size(data,2))])
        line([500 500],[-1.5 1.5],'color','r')
        line([2500 2500],[-1.5 1.5],'color','r')
        axis([0 4500 -1.5 1.5])
        
        cnt=cnt+1;
    end
end


%% Plot for observations - overlay of mean responses - black flash and white flash

col='brg'
conditions={'control','apb','wash'};

figure
contr='black';
for i=1:3 % goes through conditions
    cnt=1;
    for nd=[6 5 4]        
        subplot(2,3,cnt) % goes through ND
        [data,~]=goodCellsSelection(datarun_nd, offCells, nd, contr, conditions{i}, 'threshRev', reversed,'threshStable', stable);
        
        hold on
        plot(mean(data,2),col(i),'linewidth',2)

        title(['ND',int2str(nd),', ',contr,' flash'])
        cnt=cnt+1;
    end
end

contr='white';
for i=1:3 % goes through conditions
    cnt=4;
    for nd=[6 5 4]        
        subplot(2,3,cnt) % goes through ND
        [data,~]=goodCellsSelection(datarun_nd, offCells, nd, contr, conditions{i}, 'threshRev', reversed,'threshStable', stable);

        hold on
        plot(mean(data,2),col(i),'linewidth',2)
        title(['ND',int2str(nd),', ',contr,' flash'])
        cnt=cnt+1;
    end
end

%subplots options
for i=1:6
    subplot(2,3,i)
    legend('Contr','APB','wash')
    line([500 500],[-0.75 1.05],'color','k')
    line([2500 2500],[-0.75 1.05],'color','k')
    line([0 4500],[0 0],'color','k')
    axis([0 4500 -0.75 1.05])
end



%%
%%%% plot for single responses


contr='black';
nd=6;

cond='control';
[data,cellIDs]=goodCellsSelection(datarun_nd, offCells, nd, contr, cond, 'threshRev', reversed,'threshStable', stable);

cond='apb';
[data1,cellIDs1]=goodCellsSelection(datarun_nd, offCells, nd, contr, cond, 'threshRev', reversed,'threshStable', stable);

cond='wash';
[data2,cellIDs2]=goodCellsSelection(datarun_nd, offCells, nd, contr, cond, 'threshRev', reversed,'threshStable', stable);

totalCells=union(cellIDs,union(cellIDs1,cellIDs2));
[rows, cols]=opt_subplots(length(totalCells));

figure
for i=1:length(totalCells)
    subplot(rows,cols,i)
    hold on
    if ismember(totalCells(i),cellIDs)
        plot(data(:,cellIDs==totalCells(i)))
    end
    if ismember(totalCells(i),cellIDs1)
        plot(data1(:,cellIDs1==totalCells(i)),'r')
    end
    if ismember(totalCells(i),cellIDs2)
        plot(data2(:,cellIDs2==totalCells(i)),'g')
    end
    title(['ND',int2str(nd),', ',contr,', cell ' int2str(totalCells(i))])
    
    axis tight
    a=get(gca,'YLim');
    line([500 500],[a(1),a(2)],'color','k')
    line([2500 2500],[a(1),a(2)],'color','k')
    line([0 4500],[0,0],'color','k')
end
    

%% response quantification

responseCode=zeros(size(datarun_nd,2),3*2*3*3); % cells, response type, contrast, condition, nd

conditions={'control','apb','wash'};
contrasts={'black','white'};

coeff=2;
cnt=0;
for nd=[6 5 4] % goes through nd
    
    for j=1:3 % goes through conditions
        
        for k=1:2 % goes through contrasts
            if k==1
                rangesCorrection=0;
            else
                rangesCorrection=-2000;
            end
            
            [data,cellIDs]=goodCellsSelection(datarun_nd, offCells, nd, contrasts{k}, conditions{i}, 'threshRev', reversed,'threshStable', stable);

            for i=1:length(cellIDs) % goes through cells
                noiseSTD=std(data(1:450,i));
                if isnan(noiseSTD)
                    noiseSTD=0;
                end
                
                %ON inhibition
                if min(data((2530:2800)+rangesCorrection,i))<min(-noiseSTD*coeff,-0.07)...
                        &&min(data((2530:2800)+rangesCorrection,i))<min(data(1:450,i))*1.3
                    responseCode(cellIDs(i),1+cnt)=1;
                end
                
                %early ON excitation
                if mean(data((2600:2850)+rangesCorrection,i))>noiseSTD*coeff...
                        &&max(data((2600:2850)+rangesCorrection,i))>max(data(1:450,i))*1.3...
                        &&max(data((2600:2850)+rangesCorrection,i))>max(data((2400:2500)+rangesCorrection,i))...
                        &&(max(data((2600:2850)+rangesCorrection,i))-min(data((2550:2650)+rangesCorrection,i)))>0.07
                    responseCode(cellIDs(i),2+cnt)=1;
                end
                
                %delayed ON excitation
                if mean(data((2950:3700)+rangesCorrection,i))>noiseSTD*coeff||(max(data((2950:3700)+rangesCorrection,i))>max(data(1:450,i))*1.5&&max(data((2950:3700)+rangesCorrection,i))>0.1)
                    responseCode(cellIDs(i),3+cnt)=1;
                end
            end
            cnt=cnt+3;
        end
    end
end
  


responseCode=reshape(responseCode,111,3,2,3,3); % cells, response type, contrast(x1), condition(x2), nd(x3)


save('/Users/alexth/Desktop/old_stuff/ONresp/matrixForThomas.mat','responseCode','offCells')

%visual check

figure
% what to plot
contr='black';
nd=6;
cond='control';

%parse input
if strcmp(contr,'black')
    x1=1;
else
    x1=2;
end

if strcmp(contr,'control')
    x2=1;
elseif strcmp(contr,'apb')
    x2=2;
else
    x2=3;
end

if nd==6
    x3=1;
elseif nd==5
    x3=2;
else
    x3=3;
end

[data,cellIDs]=goodCellsSelection(datarun_nd, offCells, nd, contr, cond, 'threshRev', reversed,'threshStable', stable);

[rows, cols]=opt_subplots(length(cellIDs));
for i=1:length(cellIDs)
    subplot(rows,cols,i)
    hold on
    plot(data(:,i))

    title(int2str(responseCode(cellIDs(i),(1:3),x1,x2,x3)))
    axis tight
    a=get(gca,'YLim');
    line([500 500],[a(1),a(2)],'color','k')
    line([2500 2500],[a(1),a(2)],'color','k')
    line([0 4500],[0,0],'color','k')
end


%% OFF response black flash, control - APB

% ND6

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,2,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).black.nd6(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(500:1000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).black.nd6(:,2);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(500:1000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(1,1)=nanmean(diff(ret))
returns(1,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(1,1)=mean(diff(lat))
latencies(1,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(1,1)=mean(diff(ampl))
amplitudes(1,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))


% ND5

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,3,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).black.nd5(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(500:1000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).black.nd5(:,2);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(500:1000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(2,1)=nanmean(diff(ret))
returns(2,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))


latencies(2,1)=mean(diff(lat))
latencies(2,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(2,1)=mean(diff(ampl))
amplitudes(2,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))




% ND4

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,4,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).black.nd4(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(500:1000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+499:2500)<a,1);    
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).black.nd4(:,2);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(500:1000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 
 

returns(3,1)=nanmean(diff(ret))
returns(3,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(3,1)=mean(diff(lat));
latencies(3,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(3,1)=mean(diff(ampl))
amplitudes(3,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))

% stats
figure
subplot(1,3,1)
errorbar(returns(:,1),returns(:,2),'.')
hold on
bar(returns(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency to 1/3 of peak amplitude','control-APB'})
axis([0 4 -300 50])

subplot(1,3,2)
errorbar(latencies(:,1),latencies(:,2),'.')
hold on
bar(latencies(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency of peak amplitude','control-APB'})
axis([0 4 -300 50])

subplot(1,3,3)
errorbar(amplitudes(:,1),amplitudes(:,2),'.')
hold on
bar(amplitudes(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Peak amplitude','control-APB'})
axis([0 4 -300 50])


figure
subplot(1,2,1)
errorbar(returns(:,1),returns(:,2),'.')
hold on
bar(returns(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency to 1/3 of peak amplitude','control-APB'})
axis([0 4 -300 50])

subplot(1,2,2)
errorbar(latencies(:,1),latencies(:,2),'.')
hold on
bar(latencies(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency of peak amplitude','control-APB'})
axis([0 4 -300 50])



%% OFF response white flash, control- APB



% ND6

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,2,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).white.nd6(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(2500:3000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+2499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+2499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).white.nd6(:,2);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(2500:3000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+2499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+2499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(1,1)=nanmean(diff(ret))
returns(1,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(1,1)=mean(diff(lat))
latencies(1,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(1,1)=mean(diff(ampl))
amplitudes(1,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))


% ND5

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,3,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).white.nd5(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(2500:3000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+2499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+2499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).white.nd5(:,2);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(2500:3000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+2499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+2499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(2,1)=nanmean(diff(ret))
returns(2,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))


latencies(2,1)=mean(diff(lat))
latencies(2,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(2,1)=mean(diff(ampl))
amplitudes(2,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))




% ND4

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,4,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).white.nd4(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(2500:3000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+2499:4500)<a,1);    
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+2499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+2499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).white.nd4(:,2);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(2500:3000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+2499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+2499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 
 

returns(3,1)=nanmean(diff(ret))
returns(3,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(3,1)=mean(diff(lat))
latencies(3,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(3,1)=mean(diff(ampl))
amplitudes(3,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))




%% OFF response black flash, control - wash

% ND6

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,2,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).black.nd6(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(500:1000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).black.nd6(:,3);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(500:1000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(1,1)=nanmean(diff(ret))
returns(1,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(1,1)=mean(diff(lat))
latencies(1,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(1,1)=mean(diff(ampl))
amplitudes(1,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))


% ND5

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,3,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).black.nd5(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(500:1000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).black.nd5(:,3);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(500:1000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(2,1)=nanmean(diff(ret))
returns(2,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))


latencies(2,1)=mean(diff(lat))
latencies(2,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(2,1)=mean(diff(ampl))
amplitudes(2,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))




% ND4

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,4,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).black.nd4(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(500:1000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+499:2500)<a,1);    
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).black.nd4(:,3);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(500:1000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 
 

returns(3,1)=nanmean(diff(ret))
returns(3,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(3,1)=mean(diff(lat))
latencies(3,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(3,1)=mean(diff(ampl))
amplitudes(3,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))

% stats
figure
subplot(1,3,1)
errorbar(returns(:,1),returns(:,2),'.')
hold on
bar(returns(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency to 1/3 of peak amplitude','control-wash'})
axis([0 4 -300 50])

subplot(1,3,2)
errorbar(latencies(:,1),latencies(:,2),'.')
hold on
bar(latencies(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency of peak amplitude','control-wash'})
axis([0 4 -300 50])

subplot(1,3,3)
errorbar(amplitudes(:,1),amplitudes(:,2),'.')
hold on
bar(amplitudes(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Peak amplitude','control-wash'})
axis([0 4 -300 50])


figure
subplot(1,2,1)
errorbar(returns(:,1),returns(:,2),'.')
hold on
bar(returns(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency to 1/3 of peak amplitude','control-wash'})
axis([0 4 -300 50])

subplot(1,2,2)
errorbar(latencies(:,1),latencies(:,2),'.')
hold on
bar(latencies(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency of peak amplitude','control-wash'})
axis([0 4 -300 50])



%% OFF response white flash, control- wash



% ND6

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,2,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).white.nd6(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(2500:3000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+2499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+2499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).white.nd6(:,3);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(2500:3000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+2499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+2499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(1,1)=nanmean(diff(ret))
returns(1,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(1,1)=mean(diff(lat))
latencies(1,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(1,1)=mean(diff(ampl))
amplitudes(1,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))


% ND5

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,3,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).white.nd5(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(2500:3000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+2499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+2499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).white.nd5(:,3);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(2500:3000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+2499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+2499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 

returns(2,1)=nanmean(diff(ret))
returns(2,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))


latencies(2,1)=mean(diff(lat))
latencies(2,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(2,1)=mean(diff(ampl))
amplitudes(2,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))




% ND4

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,4,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff'
    data=datarun_nd(i).white.nd4(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(2500:3000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+2499:4500)<a,1);    
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    subplot(6,6,cnt)
    plot(data)
    hold on
    plot(lat(1,cnt)+2499,ampl(1,cnt),'k*')
    plot(ret(1,cnt)+2499,ampl(1,cnt)/3,'r*')
    
    data=datarun_nd(i).white.nd4(:,3);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(2500:3000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+2499:4500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')
    plot(lat(2,cnt)+2499,ampl(2,cnt),'k*')
    plot(ret(2,cnt)+2499,ampl(2,cnt)/3,'r*')
    
    cnt=cnt+1;
end
 
 

returns(3,1)=nanmean(diff(ret))
returns(3,2)=nanstd(diff(ret))/sqrt(sum(~isnan(diff(ret))))
[p,h]=ranksum(ret(1,:),ret(2,:))

latencies(3,1)=mean(diff(lat))
latencies(3,2)=std(diff(lat))/sqrt(size(lat,2))
[p,h]=ranksum(lat(1,:),lat(2,:))

amplitudes(3,1)=mean(diff(ampl))
amplitudes(3,2)=std(diff(ampl))/sqrt(size(ampl,2))
[p,h]=ranksum(ampl(1,:),ampl(2,:))

% stats
figure
subplot(1,3,1)
errorbar(returns(:,1),returns(:,2),'.')
hold on
bar(returns(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency to 1/3 of peak amplitude','control-wash'})
axis([0 4 -300 50])

subplot(1,3,2)
errorbar(latencies(:,1),latencies(:,2),'.')
hold on
bar(latencies(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency of peak amplitude','control-wash'})
axis([0 4 -300 50])

subplot(1,3,3)
errorbar(amplitudes(:,1),amplitudes(:,2),'.')
hold on
bar(amplitudes(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Peak amplitude','control-wash'})
axis([0 4 -300 50])


figure
subplot(1,2,1)
errorbar(returns(:,1),returns(:,2),'.')
hold on
bar(returns(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency to 1/3 of peak amplitude','control-wash'})
axis([0 4 -300 50])

subplot(1,2,2)
errorbar(latencies(:,1),latencies(:,2),'.')
hold on
bar(latencies(:,1))
set(gca,'XTick',1:3,'XTickLabel',{'ND6','ND5','ND4'})
title({'Delta latency of peak amplitude','control-wash'})
axis([0 4 -300 50])


%% Examples


% OFF response shortening (black flash, ND6)

threshRev=0.6;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,2,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));
cnt=1;
figure
for i=trustOff(20)
    data=datarun_nd(i).black.nd6(:,1);
    data=data-mean(data(1:450));
    [ampl(1,cnt) lat(1,cnt)]=max(data(500:1000));
    a=ampl(1,cnt)/3;
    a=find(data(lat(1,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(1,cnt)-1;
        ret(1,cnt)=a;
    else
        ret(1,cnt)=NaN;
    end

    plot(data)
    hold on

    
    data=datarun_nd(i).black.nd6(:,2);
    data=data-mean(data(1:450));
    [ampl(2,cnt) lat(2,cnt)]=max(data(500:1000));
    a=ampl(2,cnt)/3;
    a=find(data(lat(2,cnt)+499:2500)<a,1);
    if ~isempty(a)
        a=a+lat(2,cnt)-1;
        ret(2,cnt)=a;
    else
        ret(2,cnt)=NaN;
    end
    
    plot(data,'g')

end

legend('Control','APB')
plot(lat(1,cnt)+499,ampl(1,cnt),'k*','markersize',10)
plot(ret(1,cnt)+499,ampl(1,cnt)/3,'r*','markersize',10)

plot(lat(2,cnt)+499,ampl(2,cnt),'k*','markersize',10)
plot(ret(2,cnt)+499,ampl(2,cnt)/3,'r*','markersize',10)

line([0 4500], [0 0],'color','k')
axis tight
a=get(gca,'YLim');
line([500 500], [a(1) a(2)*1.1],'color','k')
line([2500 2500], [a(1) a(2)*1.1],'color','k')
tmp=(ampl(1,cnt)/3+ampl(2,cnt)/3)/2;
line([ret(1,cnt)+499 ret(2,cnt)+499],[tmp tmp],'color','r')
axis tight


% % checkNames
threshRev=0.3;
threshStable=0.3;

goodOff=intersect(offCells,find(sum(squeeze(ifStable(:,2,:))>threshStable,2)==3));
trustOff=intersect(goodOff,find(ifReversable(:,1)>threshRev));

numel(trustOff)

for i=trustOff'
    data=datarun_nd(i).name
end