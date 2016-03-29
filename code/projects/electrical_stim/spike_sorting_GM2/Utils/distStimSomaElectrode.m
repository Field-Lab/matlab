function [isSame distance] = distStimSomaElectrode(templates,stimElec)
%finds the distance between stimulating elecrode and the 'somatic'
%electrodes, which gives a proxy for somatic position
%Gonzalo Mena 03/16
load arrayPositions512
for g         = 1:length(templates)
    [a b]     = sort(max(abs(templates{g}')),'descend');
    isSame(g) = stimElec==b(1);
    distance(g) = norm(positions(stimElec,:)-positions(b(1),:));
end
        