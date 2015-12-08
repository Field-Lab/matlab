function plot_raster(raster, start_time, end_time)
% plot_raster(raster, start_time, end_time)
% raster: nx1 cells, one cell for a trial
%
% xyao
% 2013-12-16


for j = 1:length(raster)
    SpikeTime = raster{j};
    SpikeTime = SpikeTime';
    X = [SpikeTime; SpikeTime];
    Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
    line(X, Y, 'color', 'b');
    axis([start_time, end_time,0,length(raster)])
    hold on
end

