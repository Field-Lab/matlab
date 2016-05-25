function plot_sta(stacell, t_frame, fh)
%PLOT_STA Plots a greyscale STA
%
%  PLOT_STA(STACELL, T_FRAME) plots STACELL with a new frame every T_FRAME
%
%  PLOT_STA(..., FH) plots in figure with handle fh

if nargin == 3
    fh = -1;
end
    
[minsta, maxsta] = find_range_sta_cell(stacell);

if fh > 0
    figure(fh); 
end
colormap(gray);
caxis([minsta, maxsta]);

for k=1:length(stacell)
    imagesc(stacell{k})
    caxis([minsta, maxsta])
    pause(t_frame)
end

end % find_range_sta_cell