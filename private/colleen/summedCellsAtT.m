function [summedCellsR, summedCellsL] = summedCellsAtT(time, flt_rsp2, flt_rsp2_shiftedRight, flt_rsp2_shiftedLeft, indices1, spks_2)
            time = num2cell((repmat(time, size(indices1,2),1))) ;
            time_cell = cell(size(indices1'));
            time_cell = cellfun(@(x,t) t, time_cell, time, 'UniformOutput', false);
            emptyCells = cellfun(@isempty,spks_2);
            everyCellAtTRight = cell2mat(cellfun(@(x,t) x(t), flt_rsp2_shiftedRight,time_cell, 'UniformOutput', false));
           everyCellAtTLeft = cell2mat(cellfun(@(x,t) x(t), flt_rsp2_shiftedLeft,time_cell, 'UniformOutput', false));

               everyCellAtTLeft(emptyCells) = 0;
               everyCellAtTRight(emptyCells) = 0;               
             summedCellsR = sum(everyCellAtTRight);
             summedCellsL = sum(everyCellAtTLeft);
           