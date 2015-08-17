function [predictedResponse,lam,nr,dr] = predict_divGainASM(fitASM,X,Xtilde,interval,f,nTrials)
nSU=length(fitASM.params);
params = fitASM.params;


totalInp = zeros(1,size(X,2));
for isu=1:nSU
totalInp = totalInp +  f((params{isu}.k.value'*X + params{isu}.b.value)./(params{isu}.a.value'*Xtilde + params{isu}.sig.value) + params{isu}.c.value);
end

nr = params{isu}.k.value'*X + params{isu}.b.value;
dr = (params{isu}.a.value'*Xtilde + params{isu}.sig.value);


lam = totalInp *interval/120;


predictedResponse = zeros(nTrials,length(lam));

for itrial=1:nTrials
predictedResponse(itrial,:)=poissrnd(lam);    
end


end