function [Gibbs Log] = HdeleteSpikesBeginning(input,Gibbs,Log,maxConds)

nNeurons = Gibbs.params.nNeurons;
I        = Gibbs.params.I;
Thres    = input.params.Heuristic.ActivationThres;


for n=1:nNeurons
    contLog(n) = Log(n).params.contLogHeuristic;
    for j=1:maxConds
        
        CondsDel{n} = [];
        
        if(nansum(Gibbs.variables.spikes{n}(j,:))./I(j)>Thres)
            CondsDel{n} = [CondsDel{n} j];
        end
        
    end
    if(~isempty(CondsDel{n}))
        contLog(n)=contLog(n)+1;
        Log(n).Heuristic{contLog(n)}=['Deletion of High spiking at the begining (Condition) ' num2str(CondsDel{n})];
        Log(n).params.contLogHeuristic = contLog(n);
    end
end
Gibbs = deleteSpikesAndResample(Gibbs,CondsDel);
