codebase_path = matlab_code_path; 

dataPath = '/Volumes/Analysis/2012-09-24-3/data000/data000'; 



datarun  = load_data(dataPath);

datarun  = load_neurons(datarun);

datarun  = load_sta(datarun, 'load_sta', 'all');

datarun  = load_params(datarun);

datarun  = load_ei(datarun, cellIds,'keep_java_ei','false');



    ei = datarun.ei.eis{get_cell_indices(datarun,cellIds(1))}; 
    eiAmps = max(ei,[],2) - min(ei,[],2);


switch length(eiAmps)
    case 512
        [elec_x, elec_y] = getElectrodeCoords512();
    case 519
        [elec_x, elec_y] = getElectrodeCoords519();
    case {61,64}
        [elec_x, elec_y] = getElectrodeCoords61();
end


figure

validIDs = [];

for n = 1: length(cellIds)
    ei = datarun.ei.eis{get_cell_indices(datarun,cellIds(n))}; 
    eiAmps = max(ei,[],2) - min(ei,[],2);
  
    scatter(elec_x, elec_y, eiAmps+.1, 'filled');
    hold on;
    
    [axon_x, axon_y] = weighted_axon_poly_reg(eiAmps);
 
    plot(axon_x, axon_y);
    
    yn = input('Valid fit?');
    
    if yn == 1
        validIDs = [validIDs cellIds(n)];
    end    
    clf;
end
validIDs