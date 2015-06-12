function Output = OutputResults(input,Gibbs,Log)
%OutputResults Gives all the relevant information
%   input:  input,Gibbs,Log structure after executing spike Sorting
%   Output: Output structure that contains spikes{n}(j,i), latencies{n}(j,i), Residuals 
%           Residual{j,e}(i,t) standard deviations sigma(e,j),
%           Artifact{e}(j,t), Artifact substracted data traces
%          (ResidualArtifact{j,e}(i,t)), substracted traces (ResidualSpikes{j,e}(i,t),the logistic regression fit for each neuron
%           LogisticReg(n,j). 
%           It also contains information about neurons in Output.neuronInfo
%           (templates, neuron Id, etc), information about stimulus in Output.stimInfo (amplitudes of stimulation,
%           breakpoints in the stimulating electrodes, etc) and information
%           about data traces in the recording electrodes in
%           Output.tracesInfo (including data itselt, as
%           Output.tracesInfo.data{j,e}(i,t) and breakpoints in data traces
%           (both axonal and coming from stimulating electrode)
%           Finally, it has the Log of the Heuristics for each  neuron,
%           Output.Log(n), the parameters used in Output.params and the
%           path in Output.path;
%                      
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
Output.Log  = Log;
Output.neuronInfo = input.neuronInfo;
Output.stimInfo   = input.stimInfo;
Output.tracesInfo = input.tracesInfo;
Output.params = input.params;
Output.path = input.names.path;
end
