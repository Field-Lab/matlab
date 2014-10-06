em512 = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(501);
em519 = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1501);


%%
electrodeMap = em519;

nparts = 4;
for i = 1:nparts
    regions{i}  = electrodeMap.getRegionShape(i, nparts);
    oregions{i} = electrodeMap.getRegionShape(i, nparts, electrodeMap.getPitch());
    emparts{i}  = electrodeMap.getElectrodeMap(i, nparts);
end; clear i nparts


%%
empart = emparts{4};
num_positions = empart.getNumberOfElectrodes - 1; % Don't include the first trigger electrode
positions = zeros(num_positions,2);
for pp = 1:num_positions
    positions(pp,1) = empart.getXPosition(pp);
    positions(pp,2) = empart.getYPosition(pp);
end; clear pp num_positions

core = false(length(positions), 1);
for i = 1:length(positions)
    core(i) = empart.isCoreElectrode(i);
end

axis equal;
hold on;
plot(positions(core, 1), positions(core, 2), 'b.');
plot(positions(~core, 1), positions(~core, 2), 'r.');


%%
axis equal;
hold on;

styles = {'r.', 'bo', 'gs', 'k^'};
for i = 1:length(emparts)
    electrodeMap = emparts{i};
    
    num_positions = electrodeMap.getNumberOfElectrodes - 1; % Don't include the first trigger electrode
    positions = zeros(num_positions,2);
    for pp = 1:num_positions
        positions(pp,1) = electrodeMap.getXPosition(pp);
        positions(pp,2) = electrodeMap.getYPosition(pp);
    end; clear pp num_positions
    
    plot(positions(:,1), positions(:,2), styles{i});
end; clear i r

styles = {'r--' 'b--' 'g--' 'k--'};
for i = 1:length(oregions)
    r = oregions{i};
    plot(r.xPoints(1:r.nPoints), r.yPoints(1:r.nPoints), styles{i});
end; clear i r

styles = {'r' 'b' 'g' 'k'};
for i = 1:length(regions)
    r = regions{i};
    plot(r.xPoints(1:r.nPoints), r.yPoints(1:r.nPoints), styles{i});
end; clear i r