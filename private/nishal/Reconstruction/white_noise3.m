
%finalLen from before
%noCells from before
load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_movX.mat');
load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_2.mat');

delay=100;
noCells=length(spk_coll);

%pixX=size(mov,1)/2;
%pixY=size(mov,2)/2;

%pixCol=1:1;
%pixX
%pixY
%pixCol
newTrialStart=1:3600*3:607844
mov_recons = (mov/255-mean(mov/255));
vecLen=delay*noCells;
q=zeros(vecLen,1);
P=zeros(vecLen,vecLen);
P_log=cell(56,1);
q_log=cell(56,1);
%A=zeros(nSamples,delay*noCells+1);
%b=zeros(nSamples,1);
icnt=0;

for iTrial=22:56%1800*120
    iTrial
    startIdx=newTrialStart(iTrial)+3600+1;
    endIdx = newTrialStart(iTrial)+3600*3-2*delay;
    
    for iLen=startIdx:endIdx
    icnt=icnt+1;
        a=[];% add a constant
        % a=[1] ?? 
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
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
save('/Volumes/Analysis/nishal/recons_3.mat','q','P','filter');

% Make prediction
nTestTrials=56;
mov_pred=zeros(3600,nTestTrials);
recons_idx=trainsampleLen+1:finalLen;%1:trainsampleLen%
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
        a=[a;double(spk_coll{icell}(iLen:1:iLen+delay-1))];
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
stairs(m2);
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
