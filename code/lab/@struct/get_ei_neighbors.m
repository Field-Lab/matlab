function datarun = get_ei_neighbors(datarun, num_neighbor)

if nargin < 2, num_neighbor = 1; end

array_id = datarun.ei.java_ei.getArrayID();
electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(array_id);

numelectrodes = electrodeMap.getNumberOfElectrodes() - 1;
neighbors = cell(numelectrodes, 1);
for i = 1:numelectrodes
    n = electrodeMap.getAdjacentsTo(i, num_neighbor);
    neighbors{i} = n(2:end);
end

datarun.ei.neighbors = neighbors;