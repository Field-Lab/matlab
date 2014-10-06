
% load images and datarun
run /jacob/snle/lab/Experiments/Array/Analysis/2007-09-18-4/images/image_catalog.m



% for each match, plot the axon and EI
for mm = 1:size(match_list,1)
    axon_id = match_list(mm,1);
    cell_id = match_list(mm,2);
    
    plot_ei_scroll(datarun,cell_id,'scale',3,'axon',axons{axon_id},'figure',cell_id)
end
