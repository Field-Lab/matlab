function Types=NS_ClusterAndPrint(StimElectrodes,RecElectrodes,Movies,BadElectrodes,ReadPath,FileName,WritePath);

for i=StimElectrodes
    for j=Movies
        [Types,TypesNew,EIs,Traces]=NS_PCA_test_function(i,j,RecElectrodes,BadElectrodes,ReadPath,FileName,WritePath);
        %clear;
    end
end