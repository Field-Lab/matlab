function Traces=NS_PrintOverlaidResponsesManyAmplitudes(StimElectrodes,RecElectrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath,FigureProperties);

for i=StimElectrodes
    for j=Movies
        FileName
        [Traces]=NS_PrintOverlaidResponses(i,j,RecElectrodes,BadElectrodes,ReadPath,FileName,WritePath,FigureProperties);
        %clear;
    end
end