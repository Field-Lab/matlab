
function [predictedResponse,lam] = predictGMLM_bias_lr_opnl(fitGMLM,maskedMov,nTrials,interval)

[~,lam]=predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,interval)
lam = fitGMLM.NL_op(lam) ;

predictedResponse = zeros(nTrials,length(lam));

for itrial=1:nTrials
predictedResponse(itrial,:)=poissrnd(lam);    
end

predictedResponse= predictedResponse';
end