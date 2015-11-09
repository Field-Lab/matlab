
function [predictedResponse,lam] = predictGMLM_gamma_opnl(fitGMLM,maskedMov,nTrials,interval)

[~,lam]=predictGMLM_gamma2(fitGMLM,maskedMov,nTrials,2,interval)
lam = fitGMLM.NL_op(lam) ;

predictedResponse = zeros(nTrials,length(lam));

for itrial=1:nTrials
predictedResponse(itrial,:)=poissrnd(lam);    
end

predictedResponse= predictedResponse';
end