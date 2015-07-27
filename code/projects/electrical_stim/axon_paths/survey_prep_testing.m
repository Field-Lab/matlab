

codebase_path = matlab_code_path; 

dataPath = '/Volumes/Analysis/2015-04-09-2/data001/data001'; 
%invalidEIs = [616 694 1084 1234 1877 2046 2163 2191 2508 3076 3091 3276 3443 3647 3783 3827 3846 4025 4353 4861 4984 5091 5371 5401 5404 5416 5462 5553 5808 5887 6137 6411 6413 6738 7491];

%65 5134;

cellIds = [1231];
 
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

[~, cell_n] = size(cellIds);

 stats = [];
for n = 1: cell_n
    
    subplot(ceil(cell_n / 2), (cell_n > 1) + 1, n)
    
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

