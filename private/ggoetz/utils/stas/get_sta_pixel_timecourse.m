function tcpx = get_sta_pixel_timecourse(stacell, xcoord, ycoord)
%GET_STA_PIXEL_TIMECOURSE Plots an STA
%
%  [TCPX] = GET_STA_PIXEL_TIMECOURSE(STACELL, X_COORD, Y_COORD) plots the 
%  time course of the STA at pixels X_COORD, Y_COORD.

tcpx = zeros(length(stacell), 1);

% x, y coordinates switches in the STA cell array compared to what imagesc
% does
for kk=1:length(stacell)
    tcpx(kk) = stacell{kk}(ycoord, xcoord);
end

end % find_range_sta_cell