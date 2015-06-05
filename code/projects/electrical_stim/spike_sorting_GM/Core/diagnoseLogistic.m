function Gibbs=diagnoseLogistic(Gibbs)


J  = Gibbs.params.J;
I  = Gibbs.params.I;
lambdaLogReg = Gibbs.params.lambdaLogReg;
alphaLogReg  = Gibbs.params.alphaLogReg;
nNeurons  = Gibbs.params.nNeurons;
spikes    = Gibbs.variables.spikes;
latencies = Gibbs.variables.latencies;
Probs     = Gibbs.variables.Probs;
alpha     = Gibbs.variables.alpha;
Gibbs.diagnostics.Logisticsp1GoF = chi2cdf(Gibbs.diagnostics.LogisticDeviance,J-3,'upper');

for n=1:nNeurons

  Deviance(n) = 0;
  
    for j = 1:J
        nspikes(j) = sum(spikes{n}(j,:));
        
        sum1 = 0;
        sum2 = 0;
        if(nspikes(j)>0);
        sum1 = 2*nspikes(j)*log(nspikes(j)/(I(j)*Probs(n,j)));
        end
        if(nspikes(j)<I)
        sum2 = 2*(I(j)-nspikes(j))*log((I(j)-nspikes(j))/(I(j)*(1-Probs(n,j))));
        end

        residual(n,j) = sum1+sum2;
        Deviance(n) = Deviance(n)+residual(n,j);
    end

    error(n,:) = abs(nanmean(spikes{n}')-Probs(n,:));


    [errorSorted indexes]=sort(-error(n,:));
    condSortError(n,:)=indexes;
    [residualSorted indexes]=sort(-residual(n,:));
    condSortResidual(n,:)=indexes;

    cmaxresidual= condSortResidual(n,1);

    J2=setdiff([1:J],cmaxresidual);
    covariates=[nansum(spikes{n}(J2,:)')' I(J2)'];

    [B FitInfo] = lassoglm(J2',covariates,'binomial','link','logit','Lambda',lambdaLogReg,'Alpha',alphaLogReg);

    Gibbs.diagnostics.logisticDeviance2(n,:) = FitInfo.Deviance;
end
Gibbs.diagnostics.LogisticDeviance3 = Deviance;

Gibbs.diagnostics.condSortError    = condSortError;
Gibbs.diagnostics.condSortResidual = condSortResidual;
Gibbs.diagnostics.Logisticsp2GoF   = chi2cdf(Gibbs.diagnostics.LogisticDeviance,J-1-3,'upper');

