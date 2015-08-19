codebase_path = matlab_code_path; 

dataPath = '/Volumes/Analysis/2015-04-09-2/data001/data001'; 

%invalidEIs = [616 694 1084 1234 1877 2046 2163 2191 2508 3076 3091 3276 3443 3647 3783 3827 3846 4025 4353 4861 4984 5091 5371 5401 5404 5416 5462 5553 5808 5887 6137 6411 6413 6738 7491];

%65 5134;

cellIds = [256 481 616 1084 1234 1743 2508 3076 4353 4861 5091 5192 6526 6738 7261];

goodEIs = [256 481 1743 5192 6526 7261];

amacrine = [1084 1234 2508 5091 6738 7491];

multiaxon = [616 3076 4353 4861];
 
% Load EI. 

datarun  = load_data(dataPath);

datarun  = load_neurons(datarun);

datarun  = load_sta(datarun, 'load_sta', 'all');

datarun  = load_params(datarun);

datarun  = load_ei(datarun, cellIds,'keep_java_ei','false');






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

 
 ama = zeros(71,1);
 good = zeros(71,1);
 multi = zeros(71,1);
for n = 1: cell_n
    ei = datarun.ei.eis{get_cell_indices(datarun,cellIds(n))}; 
    stats = [];
    
    for t = 1:size(ei, 2)
        stats = [stats; sum(ei(:,t))];
    end  
    
    max_amp = max(stats);
    
    %stats = stats/max_amp;
    
    if ismember(cellIds(n), goodEIs)
        good = good + stats;
    elseif ismember(cellIds(n), amacrine)
        ama = ama + stats;
    else
        multi = multi + stats;
    end
    
     plot(good);
    hold on;
    plot(ama);
    hold on;
    %plot(multi);
    
    legend('Good EIs', 'Amacrine', 'Multiple Axons');

end

    
   
    
    
