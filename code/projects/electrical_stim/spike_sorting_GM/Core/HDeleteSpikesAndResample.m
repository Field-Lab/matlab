function [Gibbs2 Log] = HDeleteSpikesAndResample(Gibbs,input,Log)



ActivationThres       = input.params.Heuristic.ActivationThres;
VarianceFactor        = input.params.Heuristic.VarianceFactor;

E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
I           = Gibbs.params.I;
nNeurons    = Gibbs.params.nNeurons;
breakPoints = input.tracesInfo.breakPoints;
dataVecJ    = input.tracesInfo.dataVecJ;
Gibbs2      = Gibbs;
prefElectrodes = input.neuronInfo.prefElectrodes;

for n=nNeurons
    contLog(n) = Log(n).params.contLogHeuristic;
end

for n=1:nNeurons
    
    
    e = prefElectrodes{n}(1);
    breakPointsAux = [breakPoints{e} J];
    for breakRange = 1:length(breakPoints{e})+1
        if(breakRange==1)
            firstCondThisRange = 1;
        else
            firstCondThisRange = breakPointsAux(breakRange-1)+1;
        end
        
        Range =[firstCondThisRange:breakPointsAux(breakRange)];
        CondsDel = detectHighSpiking(Gibbs2,input,n,breakRange);
        CondDel  = CondsDel{n};
        
        if(~isempty(CondDel))
            Log(n).Deletion(breakRange) = 1;
            contLog(n)=contLog(n)+1;
            Log(n).Heuristic{contLog(n)}=['Going to delete spikes at breakpoint range ' num2str(breakRange) ' starting at condition ' num2str(CondDel)];
            Log(n).params.contLogHeuristic = contLog(n);
            
            CondExtrapolate = [];
            CondDel0        = CondDel;
            while(CondDel<=breakPointsAux(breakRange))
                
                Gibbs2 = deleteSpikesAndResample(Gibbs2,CondsDel);
                sigmaold=Gibbs.variables.sigma(e,CondDel);
                sigmanew=Gibbs2.variables.sigma(e,CondDel);
                if(sigmanew > VarianceFactor*sigmaold)
                    contLog(n)=contLog(n)+1;
                    Log(n).Heuristic{contLog(n)}=['Deletion at Condition ' num2str(CondDel) ' (breakpoint range ' num2str(breakRange) ') increased residual variance too much. Going back to the original spikes'];
                    Log(n).params.contLogHeuristic = contLog(n);
                    
                    Gibbs2.variables.spikes{n}(CondDel,:)            = Gibbs.variables.spikes{n}(CondDel,:);
                    Gibbs2.variables.ActionPotentials{n,CondDel}     = Gibbs.variables.ActionPotentials{n,CondDel};
                    Gibbs2.variables.ArtifactE{e}(CondDel,:)         = Gibbs.variables.ArtifactE{e}(CondDel,:);
                    Gibbs2.variables.Artifact(CondDel,(e-1)*T+1:e*T) = Gibbs.variables.Artifact(CondDel,(e-1)*T+1:e*T);
                    Gibbs2.variables.sigma(n,CondDel)                = Gibbs.variables.sigma(n,CondDel);
                    [ResSpikesE ResSpikesEJ]                         = substractActionPotentials(dataVecJ,Gibbs2.variables.ActionPotentials,E);
                    Gibbs2                    = GibbsSamplerSpikesArtifact(Gibbs2);
                    break
                end
                CondDel     = CondDel+1;
                CondsDel{n} = CondDel;
            end
            
            CondsAboveThres          = find(nansum(Gibbs2.variables.spikes{n}')./I>ActivationThres);
            CondsAboveThresThisRange = intersect(CondsAboveThres,[Range(1):CondDel-1]);
            
            if(CondDel>CondDel0)
                if(isempty(CondsAboveThresThisRange))
                    if(CondDel0 == firstCondThisRange)
                        contLog(n)=contLog(n)+1;
                        Log(n).Heuristic{contLog(n)}=['Deletion led to lack of activation in breakpoint range ' num2str(breakRange) '. The first deleted condition (' num2str(CondDel0) ') is the first condition of this range. Not going to extrapolate'];
                        Log(n).params.contLogHeuristic = contLog(n);
                        
                        continue
                    else

                        %CondExtrapolate =CondDel0+1;
                        Gibbs2 = extrapolateFromCondition(Gibbs2,input,n,CondExtrapolate,breakRange);
                        contLog(n)=contLog(n)+1;
                        %Log(n).Heuristic{contLog(n)}=['Deletion led to lack of activation in breakpoint range ' num2str(breakRange) '. Now Going to extrapolate starting at condition ' num2str(CondExtrapolate) ];
                         Log(n).Heuristic{contLog(n)}=['Deletion led to lack of activation in breakpoint range ' num2str(breakRange) '. Now Going to Gibbs sample ' num2str(CondExtrapolate) ];
                        
                        Gibbs2 = GibbsSamplerSpikesArtifact(Gibbs2);

                        Log(n).params.contLogHeuristic = contLog(n);
                        
                    end
                else
                    if(isequal(CondsAboveThresThisRange,firstCondThisRange))
                        CondExtrapolate = firstCondThisRange+1;
                        Gibbs2 = extrapolateFromCondition(Gibbs2,input,n,CondExtrapolate,breakRange);
                        contLog(n)=contLog(n)+1;
                        Log(n).Heuristic{contLog(n)}=['Deletion didnt supress all activation in breakpoint range ' num2str(breakRange) '. Now Going to extrapolate starting at condition ' num2str(CondExtrapolate) ' because the first activated condition was the first of this breakpoint range'];
                        Log(n).params.contLogHeuristic = contLog(n);
                        
                    else
                        CondExtrapolate = CondsAboveThresThisRange(end);
                        Gibbs2 = extrapolateFromCondition(Gibbs2,input,n,CondExtrapolate,breakRange);
                        contLog(n)=contLog(n)+1;
                        Log(n).Heuristic{contLog(n)}=['Deletion didnt supress all activation in breakpoint range ' num2str(breakRange) '. Now Going to extrapolate starting at condition ' num2str(CondExtrapolate) ];
                        Log(n).params.contLogHeuristic = contLog(n);
                        
                    end
                end
            else
                CondsExtrapolate    = CondsDel;
                CondsExtrapolate{n} = CondDel0;
                MaxCond        = breakPointsAux(breakRange);
                
                if(CondDel0 == firstCondThisRange)
                    CondsExtrapolate{n} = CondDel0+1;
                    contLog(n)=contLog(n)+1;
                    Log(n).Heuristic{contLog(n)}=['Deletion at breakpoint range ' num2str(breakRange) ' stopped at the first attempt (Condition ' num2str(CondDel0) ', the first of the breakpoint range). Now going to resample to look for missing spikes, starting condition ' num2str(CondDel0+1)];
                    Log(n).params.contLogHeuristic = contLog(n);
                    
                else
                    contLog(n)=contLog(n)+1;
                    Log(n).Heuristic{contLog(n)}=['Deletion at breakpoint range ' num2str(breakRange) ' stopped at the first attempt (Condition ' num2str(CondDel0) '). Now going to resample to look for missing spikes, starting condition ' num2str(CondDel0)];
                    Log(n).params.contLogHeuristic = contLog(n);
                    
                end
                Gibbs2         = HResampleIfLackOfSpiking(Gibbs2,input,Log,n,CondsExtrapolate,MaxCond);
            end
            
            
        else
            Log(n).Deletion(breakRange) = 0;
            
        end
    end
    
end



