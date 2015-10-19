load('axontrace_matches.mat')

codebase_path = matlab_code_path; 

dataPath = '/Volumes/Analysis/2007-09-18-4/data002-nwpca-duplicates/data002/data002'; 

%65 5134;

% match_list = [9 2296 ;11 2913; 14 2348; 15 2342; 17 3559; 41 1326; 44 4446; 80 413];
match_list = [11 2913; 17 3559;41 1326; 44 4446];
cellIds = match_list(:,2);
 
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


[match_n, ~] = size(match_list);

estimates = {};
 stats = [];
for n = 1: match_n
    
    %subplot(ceil(match_n/2), 2, n)
    figure
    
    
    ei = datarun.ei.eis{get_cell_indices(datarun,cellIds(n))}; 

    eiAmps = max(ei,[],2) - min(ei,[],2);

    [axon_x, axon_y, ~, ~, ~, ~, res] = weighted_axon_poly_reg(eiAmps);
 
    res;
    axon_x = axon_x';
    axon_y = axon_y';
    
    axon_xy = [axon_x axon_y];
    
    estimates{n} = axon_xy;
    
    traced = axons{match_list(n, 1)};

    scatter(elec_x, elec_y, 2*(eiAmps+.1), 'filled', 'MarkerFaceColor', 'black');
    hold on;
    plot(axon_x, axon_y, 'r');
    
    plot(traced(:,1), traced(:,2),'b');
    
    if n == 1
        legend('Electrode amplitude', 'Estimated', 'Traced');
    end    
    axis off
end

