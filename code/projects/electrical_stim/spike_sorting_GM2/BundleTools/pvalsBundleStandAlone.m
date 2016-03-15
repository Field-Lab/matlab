function [pval]=pvalsBundleStandAlone(Art,patternNo,findBundleTimes,nNeighborsExclude)
%Gonzalo Mena, 3/2016
[Res]=ResidualsElectrodeSimple(Art,patternNo,findBundleTimes);


    
    
pats=getNeighbors(patternNo,nNeighborsExclude);
nConds=size(Art,1);

for k=1:nConds
    
[h p]=kstest((log(Res(k,setdiff([1:512],pats)))-nanmean(log(Res(k,setdiff([1:512],pats)))))./nanstd((log(Res(k,setdiff([1:512],pats))))));

pval(k)=log(p);
end
    pval(1)=NaN;
   