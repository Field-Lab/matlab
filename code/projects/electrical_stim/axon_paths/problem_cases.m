codebase_path = matlab_code_path; 

%dataPath = '/Volumes/Analysis/2012-04-13-5/data001/data001';

%dataPath = '/Volumes/Analysis/2015-04-09-2/data001/data001';

dataPath = '/Volumes/Analysis/2012-09-24-3/data000/data000';

%65 5134;


%cellIds = [1821 2508 3076 4025 4353 5091 5206];
 
%cellIds = 4741;

cellIds = 3457;
% Load EI. 

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






estimates = {};
 stats = [];
for n = 1: length(cellIds)
    
    
    ei = datarun.ei.eis{get_cell_indices(datarun,cellIds(n))}; 

    eiAmps = max(ei,[],2) - min(ei,[],2);

    weighted_axon_poly_reg(eiAmps, 'plot', 1);
 
    disp(cellIds(n));
    
end