%plantain
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data000/data000/data000');

% peach
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data000/data000');

% apple
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data000/data000');

% blueberry
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data000/data000');


% load data
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta','all');
datarun = set_polarities(datarun);

cell_types = {1,2,3,4};

datarun = get_sta_summaries(datarun, cell_types);

for ct = 1:4
    
    cell_indices = get_cell_indices(datarun, {ct});
    
    num_cells = length(cell_indices);
    
    clear temp_tc_matrix
    for cc = 1:num_cells
        temp_tc = datarun.stas.time_courses{cell_indices(cc)};
        summed_tc = sum(temp_tc,2);
        
        temp_tc_matrix(:,cc) = summed_tc;
    end

    tc_structure(ct).tcs = temp_tc_matrix;
end
        


save plantain_tcs tc_structure

