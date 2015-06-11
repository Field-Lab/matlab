function Output = OutputResults(input,Gibbs,Log)
%OutputResults Gives all the relevant information
%   input: input,Gibbs,Log structure after executing spike Sorting
%   Output: Output structure that contains original data traces Data{j,e}(i,t),
%   spikes{n}(j,i), latencies{n}(j,i), Residuals Residual{j,e}(i,t) standard deviations sigma(e,j),
%   Artifact{e}(j,t), Artifact substracted data traces (ResidualArtifact{j,e}(i,t)) and spike
%   substracted traces (ResidualSpikes{j,e}(i,t),the logistic regression fit for each neuron
%   LogisticReg(n,j). Finally, it has the Log of the Heuristics for each
%   neuron, Log(n)
%   
Output.spikes = Gibbs.variables.spikes;
Output.latencies = Gibbs.variables.latencies;
Output.sigma  = Gibbs.variables.sigma;
Output.Artifact = Gibbs.variables.ArtifactE;
Output.Residual = Gibbs.variables.Residuals;

E = input.tracesInfo.E;
J = input.tracesInfo.J;
I = input.tracesInfo.I;

for e=1:E
    for j=1:J
Output.ResidualSpikes{j,e} = Gibbs.variables.ResSpikesEJ{e}{j};
    end
end

for e=1:E
    for j=1:J
        Output.ResidualArtifact{j,e} = input.tracesInfo.data{j,e}-repmat(Gibbs.variables.ArtifactE{e}(j,:),I(j),1);
    end
end

Output.LogisticReg = Gibbs.variables.Probs;
Output.Data = input.tracesInfo.data;
Output.Log  = Log;
end
