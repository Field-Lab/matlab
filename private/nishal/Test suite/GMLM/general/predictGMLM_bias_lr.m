
function [predictedResponse,lam,kx] = predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,interval)


filters = fitGMLM.Linear.filter;
bias = fitGMLM.Linear.bias;
nFrontEnds = length(filters);
mu=fitGMLM.mu;

mov_filtered=maskedMov;

kx=cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered + bias{ifilter});
end


binFrame=1;
lam=kx{1};
for ik=2:length(kx)
lam=lam+kx{ik};
end
lam=lam+mu;
l2=repmat(lam,[binFrame,1]);
lam=l2(:);
lam=lam*(interval/(120*binFrame)); 

predictedResponse = zeros(nTrials,length(lam));

for itrial=1:nTrials
predictedResponse(itrial,:)=poissrnd(lam);    
end

predictedResponse= predictedResponse';
end