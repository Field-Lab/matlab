function neighbors=get_ei_neighbors(electrode,array,nr_neighbor)

if nargin<2
    ar=512;
end

if nargin<3
    nr_neighbor=1;
end

if array==61
    ar=1;
end
if array==512
    ar=505;
end
if array==519
    ar=1600;
end


electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ar);  
neighbors = electrodeMap.getAdjacentsTo(electrode, nr_neighbor);
