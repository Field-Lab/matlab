
function [predictedResponse,lam] = predictGMLM_bias_lr_opnl(fitGMLM,maskedMov,nTrials,interval)
dt = 1/120;
[~,lam]=predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,interval);
lam = fitGMLM.NL(lam/dt)*dt ;
%lam(lam>0.5)=0.5;

predictedResponse = zeros(nTrials,length(lam));

for itrial=1:nTrials
predictedResponse(itrial,:)=poissrnd(lam);    
end

predictedResponse= predictedResponse';
end