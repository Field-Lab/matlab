
%finalLen from before
%noCells from before
%load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_movX.mat');
load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_2.mat');
newTrialStart=1:3600*3:607844;


delay=120;
noCells=length(spk_coll);

% Do smoothing 
smoothParam=40;
for icell=1:noCells
spk_coll{icell}=conv(spk_coll{icell},ones(smoothParam,1)/smoothParam,'same');
end
mov=conv(mov,ones(smoothParam,1)/smoothParam,'same');

%pixX=size(mov,1)/2;
%pixY=size(mov,2)/2;

%pixCol=1:1;
pixX
pixY
%pixCol


mov_recons = (mov/255-mean(mov/255));
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
    startIdx=newTrialStart(iTrial)+3600+1;
    endIdx = newTrialStart(iTrial)+3600*3-2*delay;
    
    for iLen=startIdx:step_sz:endIdx
    icnt=icnt+1;
        a=[];% add a constant
        % a=[1] ?? 
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen:step_sz:iLen+step_sz*delay-1))];
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
        a=[a;double(spk_coll{icell}(iLen:step_sz:iLen+step_sz*delay-1))];
        end
        a=full(a);
    
        mov_pred(icnt,iTrial)= filter_use'*a;

    end
end

m1= mov_recons(1:1:3600);
m2= mov_pred;
figure;
stairs(m1,'b');
hold on
stairs(m2(:,1:3));
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
icell=1;
load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_2.mat');
trialData=zeros(3600,nTrials);

for iTrial=1:nTrials
indx=newTrialStart(iTrial)+(iTrial-1):newTrialStart(iTrial)+3600-1+(iTrial-1);
trialData(:,iTrial)=spk_coll{icell}(indx);

end

figure;
plotSpikeRaster(trialData'>0,'PlotType','vertline');

