function [Gibbs Log] = HaddSpikesResample(Gibbs,input,Log)



E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
I           = Gibbs.params.I;
nNeurons    = Gibbs.params.nNeurons;
breakPoints = input.tracesInfo.breakPoints;
prefElectrodes = input.neuronInfo.prefElectrodes;

LowThres   = input.params.Heuristic.LowThres;
VeryLowThres   = input.params.Heuristic.VeryLowThres;
HighThres      = input.params.Heuristic.HighThres;
ActivationThres = input.params.Heuristic.ActivationThres;

for n=1:nNeurons
    contLog(n) = Log(n).params.contLogHeuristic;
end



for n=1:nNeurons
    
    e = prefElectrodes{n}(1);
    breakPointsAux = [breakPoints{e} J];
    if(length(breakPointsAux) == 1)
        continue
    else
        for breakRange = 2:length(breakPointsAux)
            lastCondPrevRange  = breakPointsAux(breakRange-1);
            firstCondThisRange = breakPointsAux(breakRange-1)+1;
            Range              = [firstCondThisRange:breakPointsAux(breakRange)];
            rangeLength = breakPointsAux(breakRange)-breakPointsAux(breakRange-1);
            if(rangeLength ==1 && nansum(Gibbs.variables.spikes{n}(firstCondThisRange,:))/I(firstCondThisRange)>HighThres)
                for n2=1:nNeurons
                    CondsDel{n2}=[];
                end
                CondsDel{n}=firstCondThisRange;
                Gibbs = deleteSpikesAndResample(Gibbs,CondsDel);
                
            end
            spikes    = Gibbs.variables.spikes;
            latencies = Gibbs.variables.latencies;
            
            spikeProb = nansum(spikes{n}(firstCondThisRange,:))/nansum(I(firstCondThisRange));

            if(spikeProb<LowThres&&nansum(spikes{n}(lastCondPrevRange,:))/I(:,lastCondPrevRange)>HighThres)

           
                contLog(n)=contLog(n)+1;
                Log(n).Heuristic{contLog(n)}=['Lack of activity in breakpoint range ' num2str(breakRange) ' follows activation in last condition of previous range. Going to add spikes in all trials from condition ' num2str(Range(1)) ' to condition ' num2str(Range(end))];
                Log(n).params.contLogHeuristic = contLog(n);
                latposPrev = latencies{n}(lastCondPrevRange,:);
                latposThis = latencies{n}(Range,:);
                latposThis = latposThis(:);
                latpos     = [latposPrev latposThis(:)'];
                latpos     = latpos(latpos>0);
                medianLatency = floor(nanmedian(latpos));
                for r=Range
                    Gibbs.variables.latencies{n}(r,1:I(r)) = medianLatency;
                    Gibbs.variables.spikes{n}(r,1:I(r))    = 1;
                end
                Gibbs = UpdateVariablesChangeSpikes(Gibbs,Gibbs.variables.latencies,n);
                Gibbs = GibbsSamplerSpikesArtifact(Gibbs);
                
            end
            
            
        end
    end
end