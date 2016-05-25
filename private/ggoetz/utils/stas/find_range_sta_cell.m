function [minval, maxval] = find_range_sta_cell(stacell)
%FIND_RANGE_STA_CELL Computes the min and max values of an sta cell
%
%  [MINVAL, MAXVAL] = FIND_RANGE_STA_CELL(STACELL) will return the min and
%  max of the STA STACELL, stored in an array of cells.

minval = stacell{1}(1);
maxval = stacell{1}(1);
for kk=1:length(stacell)
    cmin = min(stacell{kk}(:));
    cmax = max(stacell{kk}(:));
    
    if cmin < minval
        minval = cmin;
    end
    if cmax > maxval
        maxval = cmax;
    end
end

end % find_range_sta_cell