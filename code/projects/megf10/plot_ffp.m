function [] = plot_ffp(raster_ff, raster_ff_all, cell_idx, duration)

d = duration;
XX = 0.05:0.1:4*d-0.05;
% a = [0 d/2 d/2 3*d/2 3*d/2 5*d/2 5*d/2 7*d/2 7*d/2 4*d]';
% b = [1 1 2 2 1 1 0 0 1 1]';

a = [0 d d 2*d 2*d 3*d 3*d 4*d]';
b = [1 1 2 2 1 1 0 0]';
idx = find(~cellfun(@isempty,raster_ff),1); % get the idx of 1st non-empty cell in raster
trialn = length(raster_ff{idx});
for j = 1:length(raster_ff{cell_idx})
    SpikeTime = raster_ff{cell_idx}{j};
    SpikeTime = SpikeTime';
    X = [SpikeTime; SpikeTime];
    Y = [ones(1, length(SpikeTime))*(j-0.9); ones(1, length(SpikeTime))*j];
    line(X, Y, 'color', 'b');
    xlim([0 4*d])
    hold on
end   
line(a, b+trialn+2)
h = hist(raster_ff_all{cell_idx}, XX);
h = h/max(h)*3;
plot(XX, h+trialn+7)
xlabel('time/sec')
end