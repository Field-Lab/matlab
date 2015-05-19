
function [predictedResponse,lam] = predictGMLM_gamma2(fitGMLM,maskedMov,nTrials,gamma,interval)
filters = fitGMLM.Linear.filter;

nFrontEnds = length(filters);
mu=fitGMLM.mu;

mov_filtered=maskedMov;

kx=cell(nFrontEnds,1);
kkx = cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kkx{ifilter}= filters{ifilter}'*mov_filtered;
    kx{ifilter} = (kkx{ifilter}.*(kkx{ifilter}>0)).^gamma;
end


sub_bin=10;
lam=kx{1};
for ik=2:length(kx)
lam=lam+kx{ik};
end
lam=lam+mu;
l2=repmat(lam,[sub_bin,1]);
lam=l2(:);
lam=lam*(interval/(120*sub_bin));

predictedResponse = zeros(nTrials,length(lam));

for itrial=1:nTrials
predictedResponse(itrial,:)=poissrnd(lam);    
end

predictedResponse= predictedResponse';
end