function [isSame distance distNeigh] = distStimSomaElectrode(templates,neuronIndexs,stimElec)
%finds the distance between stimulating elecrode and the 'somatic'
%electrodes, which gives a proxy for somatic position
%Gonzalo Mena 03/16
load arrayPositions512
for g         = 1:length(neuronIndexs)
    index=neuronIndexs(g);
    [a b]     = sort(max(abs(templates{index}')),'descend');
    isSame(g) = stimElec==b(1);
    distance(g) = norm(positions(stimElec,:)-positions(b(1),:));
    cont=0;
    while(true)
        els=getNeighbors(stimElec,cont);
        if(~isempty(find(els==b(1))))
            distNeigh=cont;
            return
        else
            cont=cont+1;
        end
    end
end
        