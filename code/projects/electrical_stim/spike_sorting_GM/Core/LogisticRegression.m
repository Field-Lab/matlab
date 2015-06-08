function Gibbs = LogisticRegression(Gibbs)


J  = Gibbs.params.J;
I  = Gibbs.params.I;
lambdaLogReg = Gibbs.params.lambdaLogReg;
alphaLogReg  = Gibbs.params.alphaLogReg;
nNeurons  = Gibbs.params.nNeurons;
spikes    = Gibbs.variables.spikes;
latencies = Gibbs.variables.latencies;


for n=1:nNeurons

    [B FitInfo] = lassoglm([[1:J]'],[nansum(spikes{n}')' I'],'binomial','link','logit','Lambda',lambdaLogReg,'Alpha',alphaLogReg);
    alpha(n,:)  = [FitInfo.Intercept B];
    Deviance    = FitInfo.Deviance;
    linear      = alpha(n,1)+alpha(n,2)*[1:J];
    probs(n,:)  = 1./(1+exp(-linear));
    if(alpha(n,2)<0)
        probs(n,1:J) = nansum(nansum(spikes{n}))/sum(I);
    end
    
    Gibbs.diagnostics.LogisticDeviance(n,:) = Deviance;
    
end

Gibbs.variables.alpha = alpha;
Gibbs.variables.Probs = probs;

Gibbs = diagnoseLogistic(Gibbs);

