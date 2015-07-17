%load('axontrace_matches.mat')

codebase_path = matlab_code_path; 

dataPath = '/Volumes/Analysis/2015-05-27-0/data000/data000'; 

%65 5134;

cellIds = [2630 2223 2717];
 
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


figure

[~, cell_n] = size(cellIds)

 stats = [];
for n = 1: cell_n
    
    subplot(ceil(cell_n), 2, n)
    
    ei = datarun.ei.eis{get_cell_indices(datarun,cellIds(n))}; 

    eiAmps = max(ei,[],2) - min(ei,[],2);

    [axon_x, axon_y] = weighted_axon_poly_reg(eiAmps, 'ei_thresh', 4);
 

    scatter(elec_x, elec_y, eiAmps+.1, 'filled', 'MarkerFaceColor', 'black');
    hold on;
    plot(axon_x, axon_y, 'r');
    
    
    if n == 1
        legend('Electrode amplitude', 'Estimated', 'Traced');
    end    
    
end

