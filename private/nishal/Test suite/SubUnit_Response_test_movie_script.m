
movieLen=size(mov2,3)-(Filtlen-1);
cell_resp=zeros(movieLen,nSubunits);
sz=max(size(mov2,3)-size(subunits{1},4) + 1, 0);
for isubunit=1:nSubunits
    isubunit
cell_resp(:,isubunit)=reshape(convn(mov2,subunits{isubunit}(end:-1:1,end:-1:1,:),'valid'),[sz,1]);
end

figure;
subplot(5,1,1);
plot(cell_resp)
title('Subunit Inputs')

% 
totalInput=0*cell_resp(:,1);
for isubunit=1:nSubunits
totalInput=totalInput+subunitWeights(isubunit)*f(cell_resp(:,isubunit));
end

totalOutput=N(totalInput);
subplot(5,1,2);
plot(totalOutput);
title('Output of the model');

dt=1/120; 
%rateScale=2.5;
rateScale=10;
firingRate=rateScale*dt*totalOutput;

binnedResponses=zeros(length(firingRate),nTrials);
for iTrial=1:nTrials
binnedResponses(:,iTrial) = poissrnd(firingRate);
end

subplot(5,1,3);
plot(binnedResponses(:,1));
title('Binned Number of spikes for a trial');

avgSpkRate = sum(binnedResponses(:))/(nTrials*movieLen/120);

% Make PSTH of response!! 
addpath(genpath('../NSEM/'));
binSize=10;
[timeLog,psth_resp]=psth_calc(binnedResponses',binSize,'nonoverlap');
subplot(5,1,4);
plot(timeLog,psth_resp);
title('Response PSTH')

% Make Raster
addpath(genpath('../plotSpikeRaster_v1/'));
subplot(5,1,5)
plotSpikeRaster(binnedResponses'>0,'PlotType','vertline');


title(sprintf('Raster, Avg spike rate: %0.02f spks/sec',avgSpkRate));


figure;
subplot(3,1,1);
[E,C]=hist(cell_resp(:,1),100);
[ax,h1,h2]=plotyy(C,E,C,f(C));
% xlim(ax(1),[min(C),max(C)]);
% xlim(ax(2),[min(C),max(C)]);
xlim(ax(1),[-5,5]);
xlim(ax(2),[-5,5]);
title('Input to a sub-unit and its non-linearity');

subplot(3,1,2);
[E,C]=hist(totalInput,100);
[ax,h1,h2]=plotyy(C,E,C,N(C));
xlim(ax(1),[0,12]);
xlim(ax(2),[0,12]);
% xlim(ax(1),[min(C),max(C)]);
% xlim(ax(2),[min(C),max(C)]);
title('Total Input to the second layer');

subplot(3,1,3);
hist(totalOutput,30);
xlim([0,30])
title(sprintf('Output of the second layer neuron,Avg spike rate: %0.02f spks/sec',avgSpkRate));


figure;
for isub=1:nSubunits
subplot(nSubunits,1,isub);
[E,C]=hist(cell_resp(:,isub),100);
[ax,h1,h2]=plotyy(C,E,C,f(C));
% xlim(ax(1),[min(cell_resp(:)),max(cell_resp(:))]);
% xlim(ax(2),[min(cell_resp(:)),max(cell_resp(:))]);
xlim(ax(1),[-5,5]);
xlim(ax(2),[-5,5]);
title(sprintf('Input to sub-unit %d and its non-linearity',isub));
end

% subunit input time courses
cell_resp_len=size(cell_resp,1);

figure;
subplot(2,1,1);
col='rbkm';
for isub=1:nSubunits
plot(cell_resp(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),isub),col(isub));
hold on;
end
title('Input to sub-units');

subplot(2,1,2);
plot(max(cell_resp(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:)'),'b');
hold on
plot(min(cell_resp(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:)'),'r');
hold on
plot(max(abs(cell_resp(uint16(cell_resp_len/2-100):uint16(cell_resp_len/2+100),:))'),'k');
legend('max','min','absolute max');
title('Max and min inputs to sub-units');

