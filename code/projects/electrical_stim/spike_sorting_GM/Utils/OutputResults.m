function Output = OutputResults(input,Gibbs,Log)
%   OutputResults Gives all the relevant information summary information After Authomatic Spike Sorting
%   input:  input,Gibbs,Log structures after executing spike Sorting
%   Output: Output structure that contains:
%           -Spikes:                                      Output.spikes{n}(j,i)
%           -Latencies:                                   Output.latencies{n}(j,i)
%           -Artifact estimates                           Output.Artifact{e}(j,t)
%           -Residuals: (data minus artifact and spikes)  Output.Residual{j,e}(i,t) 
%           -Standard deviations of residuals             Output.sigma(e,j),
%           -Artifact substracted data traces             Output.ResidualArtifact{j,e}(i,t))
%           -Spike substracted data traces                Output.ResidualSpikes{j,e}(i,t)
%           -Logistic regression fits for each neuron     Output.LogisticReg(n,j)
%           -Neuron Information                           Output.neuronInfo  (templates, neuron Ids, etc)
%           -Stimulus information                         Output.stimInfo  (pattern, amplitudes of stimulation, stimulating electrodes,etc)
%           -Traces information                           Output.tracesInfo(including data itselt, in Output.tracesInfo.data{j,e}(i,t), 
%                                                         recording
%                                                         electrodes, breakpoints in data traces, both axonal and coming from stimulating electrode.
%           -Log Information                              Output.Log(n)  Log of the Heuristics for each  neuron
%           -Algorithm Parameter information              Output.params all parameters that were set (by default or not) in the input structure
%           -Path                                         Output.path the pattern directory where movies where taken
%                      
%   Gonzalo Mena 06/15
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
