function [onset onsetC]=findBundleFrompValStandAlone(pvals,listAmps,detectThreshold)
%Gonzalo Mena, 3/2016
    aux=find(pvals>detectThreshold);
    comp=lastConnectedComponent(aux);
    nConds=length(pvals);
    if(isempty(aux));
    onset=listAmps(2);
    onsetC=2;
    return
    else
    if(length(comp)==1&&comp(end)==nConds)
        if(aux(1)==2)
        onset=NaN;
        onsetC=NaN;
        return
        else
            condi=2;
        end
    elseif(length(comp)==1&&comp(end)<nConds);
        condi=comp(end)+1;
    else
    if(comp(end)==nConds)
        condi=comp(end-1)+1;
    else
        condi=comp(end)+1;
    end
    end
    
    onset=listAmps(condi);
    onsetC=condi;
    end


