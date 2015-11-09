function g = compute_op_NL(fitGMLM,maskedMovdd,spksGen)

 nTrials=1;interval=1;
 %[resp,lam] = predictGMLM_bias_lr(fitGMLM,maskedMovdd,nTrials,interval);
  [resp,lam]= predictGMLM_gamma2(fitGMLM,maskedMovdd,nTrials,2,interval);
 resp=resp';lam=lam';
%  
% figure;
% scatter(lam,spksGen_hr,1*ones(length(spksGen_hr),1),'filled');hold on;

th_list=[];
for iprc=0:0.5:100
 thr = prctile(lam,iprc);
th_list= [th_list;thr]; 
end

meanl=[];meanR=[];meanPredR = [];
for ith=1:length(th_list)-1
iidx = (lam>th_list(ith)) & (lam<=th_list(ith+1));
meanl = [meanl;mean(lam(iidx))];
meanR = [meanR;mean(spksGen(iidx))];
meanPredR = [meanPredR ; mean(resp(iidx))];
end

%[fit,p]= fit_logistic(meanl,meanR);
%g = @(x) (p(2) ./ (1+exp(-p(3).*(x-p(1)))));
% 
 f=fit(meanl(~isnan(meanR)),meanR(~isnan(meanR)),'poly1');
% a=f.p1;
 g = @(x) ((f.p1)*x + (f.p2));

figure;
plot(meanl,meanR);
hold on;
plot(meanl,g(meanl));
hold on;
plot(meanl,meanPredR);
legend('Estimated NL','Fit NL','Linear case');

end