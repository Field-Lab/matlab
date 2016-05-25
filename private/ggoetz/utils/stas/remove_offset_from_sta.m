function [stacell] = remove_offset_from_sta(stacell, offset)
%REMOVE_OFFSET_FROM_STA Removes an offset from the STA
%
%  [STACELL] = FIND_RANGE_STA_CELL(STACELL, OFFSET) will remove OFFSET from
%  each of the frames of the sta cell STACELL.

for kk=1:length(stacell)
    stacell{kk} = stacell{kk} - offset;
end

end % find_range_sta_cell