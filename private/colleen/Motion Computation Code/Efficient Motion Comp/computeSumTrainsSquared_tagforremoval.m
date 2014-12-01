%This function computes the equation sum_t(sum_i(r_i(t-\Delta t))^2
function computeSumTrainsSquared(flt_rsp1, flt_rsp2, flt_rsp1_shifted, flt_rsp2_shifted, indices1);

    
 time_cell = cell(size(indices1'));
 time_cell = cellfun(@(x) 0.5, time_cell, 'UniformOutput', false);


    everyCellAtT = cell2mat(cellfun(@(x,t) x(t), flt_rsp2_shifted,time_cell, 'UniformOutput', false));
    
    summedCellsAtT = sum(everyCellAtT);