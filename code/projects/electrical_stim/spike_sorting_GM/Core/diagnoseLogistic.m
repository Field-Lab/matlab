function Gibbs=diagnoseLogistic(Gibbs)


J  = Gibbs.params.J;
I  = Gibbs.params.I;
lambdaLogReg  = Gibbs.params.lambdaLogReg;
alphaLogReg   = Gibbs.params.alphaLogReg;
thresLogistic = Gibbs.params.thresLogistic;

nNeurons  = Gibbs.params.nNeurons;
spikes    = Gibbs.variables.spikes;
latencies = Gibbs.variables.latencies;
Probs     = Gibbs.variables.Probs;
alpha     = Gibbs.variables.alpha;

Gibbs.diagnostics.Logisticsp1GoF = chi2cdf(Gibbs.diagnostics.LogisticDeviance,J-3,'upper');

for n=1:nNeurons

  
    spikesaux{n} =  spikes{n};
    
  
    for j = 1:J
        
        
        nspikes(j) = sum(spikes{n}(j,:));
        
        if(nspikes(j)/I(j)<=thresLogistic(j))
            spikesaux{n}(j,1:I(j))=0;
        end
        if(nspikes(j)/I(j)>=1-thresLogistic(j))
            spikesaux{n}(j,1:I(j))=1;
        end
        
    end
    
    [B FitInfo] = lassoglm([[1:J]'],[nansum(spikesaux{n}')' I'],'binomial','link','logit','Lambda',lambdaLogReg,'Alpha',alphaLogReg);

    alpha(n,:)  = [FitInfo.Intercept B];
    Deviance    = FitInfo.Deviance;
    linear      = alpha(n,1)+alpha(n,2)*[1:J];
    probs(n,:)  = 1./(1+exp(-linear));
    if(alpha(n,2)<0)
        probs(n,1:J) = nansum(nansum(spikesaux{n}))/sum(I);
    end
    
    error(n,:) = abs(nanmean(spikesaux{n}')-Probs(n,:));
    Gibbs.diagnostics.LogisticDeviance2(n) = Deviance;
    
    [errorSorted indexes]=sort(-error(n,:));
    condSortError(n,:)=indexes;
end
  

    Gibbs.diagnostics.Logisticsp2GoF   = chi2cdf(Gibbs.diagnostics.LogisticDeviance2,J-3,'upper');
    

Gibbs.diagnostics.condSortError    = condSortError;

Gibbs.variables.ProbRobust       = probs;
