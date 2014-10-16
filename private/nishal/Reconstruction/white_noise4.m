
%finalLen from before
%noCells from before
%load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_movX.mat');
load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_2.mat');
%load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_all_data.mat');
load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_movXX.mat');
newTrialStart=1:3600*3:607844;


delay=120;
noCells=length(spk_coll);
spk_coll_smooth=cell(noCells,1);

% Do smoothing 
smoothParam=10;
for icell=1:noCells
spk_coll_smooth{icell}=conv(spk_coll{icell},ones(smoothParam,1)/smoothParam,'same');
end
mov_smooth=conv(mov,ones(smoothParam,1)/smoothParam,'same');
% ERROR WITH SMOOTHING .. BECUASE FRAMES ACCROSS THE MOVIES GET MERGED!!
% ..SHOULD NOT AFFECT MUCH THOUGH!
% BAD


%pixX=size(mov,1)/2;
%pixY=size(mov,2)/2;

%pixCol=1:1;
%pixX
%pixY
%pixCol


mov_recons = (mov_smooth/255-mean(mov_smooth/255));
vecLen=delay*noCells;
q=zeros(vecLen,1);
P=zeros(vecLen,vecLen);
nTrials=55;
P_log=cell(nTrials,1);
q_log=cell(nTrials,1);
%A=zeros(nSamples,delay*noCells+1);
%b=zeros(nSamples,1);
icnt=0;
step_sz=smoothParam;
for iTrial=1:55%1800*120
    iTrial
    startIdx=newTrialStart(iTrial)+3600+1+smoothParam;
    endIdx = newTrialStart(iTrial)+3600*3-2*delay-smoothParam;
    
    for iLen=startIdx:step_sz:endIdx
    icnt=icnt+1;
        a=[];% add a constant
        % a=[1] ?? 
    for icell=1:noCells
    a=[a;double(spk_coll_smooth{icell}(iLen:step_sz:iLen+step_sz*delay-1))];
    end
    a=full(a);
    y=mov_recons(iLen);

  %  A(iLen,:)=a';
  %  b(iLen)=y;
    
    q=q+y*a;
    P=P+a*a';
    end
    P_log{iTrial}=P/icnt;
    q_log{iTrial}=q/icnt;
end

 q=q;
 P=P;
 
filter=P\q;   
%filter2=(A\(A'\(A'*b)));

filter_use=filter;

% Make prediction
nTestTrials=62;
mov_pred=zeros(3600,nTestTrials);


for iTrial=1:nTestTrials
    iTrial
    startIdx=newTrialStart(iTrial);
    endIdx = newTrialStart(iTrial)+3600-delay;
    recons_idx=startIdx:endIdx;
    icnt=0;
    for iLen=recons_idx %1800*120
        icnt=icnt+1;
        %a=[1];
        a=[];
        for icell=1:noCells
        a=[a;double(spk_coll_smooth{icell}(iLen:step_sz:iLen+step_sz*delay-1))];
        end
        a=full(a);
    
        mov_pred(icnt,iTrial)= filter_use'*a;

    end
end

m1= mov_recons(1:3600);
m2= mov_pred;
figure;
stairs(m1,'b');
hold on
stairs(m2(:,1),'r');
%hold on;
ylim([-1,1]);
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
%xlim([1,100])
% % correlation 
% 
% %filter_bank{pixX,pixY,pixCol}.filter=filter_use;
% %filter_bank{pixX,pixY,pixCol}.correlation=corr(mov_pred(recons_idx),mov_recons(binFrames(recons_idx)));
% corr(mov_pred(recons_idx(1)+100:recons_idx(1)+2000),mov_recons(recons_idx(1)+100:recons_idx(1)+2000))
% corr(m1,m2)

%% See filter structure
filt_mat=zeros(delay,noCells);
for icell=1:noCells
filt_mat(:,icell)=filter_use((icell-1)*delay+1:icell*delay);
end

figure;
plot(sum(filt_mat.^2,1))

%% See data
icell=24;
nTrials=50
%load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_all_data.mat');
trialData_smooth=zeros(3600,nTrials);
trialData=zeros(3600,nTrials);
mov_data=zeros(3600,nTrials);
newTrialStart=1:3600*3:607844;
mov_smooth_data=zeros(3600,nTrials);

for iTrial=1:nTrials
indx=newTrialStart(iTrial):newTrialStart(iTrial)+3600-1;
trialData_smooth(:,iTrial)=spk_coll_smooth{icell}(indx);
trialData(:,iTrial)=spk_coll{icell}(indx);
mov_data(:,iTrial)=mov(indx);
mov_smooth_data(:,iTrial)=mov_smooth(indx);

end

figure;
subplot(3,1,1)
plotSpikeRaster(trialData_smooth'>0,'PlotType','vertline');
title(sprintf('Smooth , Param:%d',smoothParam));
subplot(3,1,2)
plotSpikeRaster(trialData'>0,'PlotType','vertline');
title('Original')
subplot(3,1,3)
plot(mov_data,'LineWidth',2)
hold on
plot(mov_smooth_data,'r')
title('Blue: Original Red: Smoothened');
xlim([0,3600])
