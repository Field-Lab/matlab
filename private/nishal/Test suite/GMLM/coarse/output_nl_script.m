 %% compute output non-linearity to ASM.
 
 figure;
 for nSU=1:5
     fitGMLM = fitGMLM_log{nSU};
 nTrials=1;
 [resp,lam] = predictGMLM_full(fitGMLM,maskedMovdd,nTrials);
 resp=resp';lam=lam';lam=lam/1200;
%  
% figure;
% scatter(lam,spksGen_hr,1*ones(length(spksGen_hr),1),'filled');hold on;

th_list=[];
for iprc=0:1:100
 thr = prctile(lam,iprc);
th_list= [th_list;thr]; 
end

meanl=[];meanR=[];
for ith=1:length(th_list)-1
iidx = (lam>th_list(ith)) & (lam<=th_list(ith+1));
meanl = [meanl;mean(lam(iidx))];
meanR = [meanR;mean(spksGen_hr(iidx))];
end

plot(meanl,meanR);
hold on;
 end


[XX,NN] =ecdf(lam);
hold on;
plotyy(meanl,meanR,NN(NN<=max(meanl)),XX(NN<=max(meanl)));

