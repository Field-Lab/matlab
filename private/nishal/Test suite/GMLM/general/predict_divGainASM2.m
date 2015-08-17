function [predictedResponse,lam,nr,dr] = predict_divGainASM2(fitASM,X,Xtilde,interval,f,nTrials)
nSU=length(fitASM.params);
params = fitASM.params;


totalInp = zeros(1,size(X,2));
nr = zeros(length(totalInp),nSU);
dr = zeros(length(totalInp),nSU);
for isu=1:nSU
totalInp = totalInp +  f((params{isu}.k.value'*X + params{isu}.b.value))./(params{isu}.a.value'*Xtilde + params{isu}.sig.value) + params{isu}.c.value;

nr(:,isu) = params{isu}.k.value'*X + params{isu}.b.value;
dr(:,isu) = (params{isu}.a.value'*Xtilde + params{isu}.sig.value);

end


lam = totalInp *interval/120;


predictedResponse = zeros(nTrials,length(lam));

for itrial=1:nTrials
predictedResponse(itrial,:)=poissrnd(lam);    
end


end