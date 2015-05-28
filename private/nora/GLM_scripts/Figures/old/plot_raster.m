function plot_raster(xval, fittedGLM, labels)

dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];

rec_rast = xval.rasters.recorded(:,1:bins);
sim_rast = xval.rasters.glm_sim(:,1:bins);

trials   = size(rec_rast,1);

hFig1=figure;
set(hFig1, 'Position', [100 100 500 250])
hold on
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'k.')
end
if labels
    xlabel('Time (seconds)')
    ylabel('Trials')
else
    axis off
end

hFig2=figure;
set(hFig2, 'Position', [100 100 500 250])
hold on
for i_trial = 1:trials
    sim1 = time(find(sim_rast(i_trial,:)));
    if length(sim1) < 4*length(rec1)
        plot(sim1, i_trial, 'r.')
    end
end
if labels
    xlabel('Time (seconds)')
    ylabel('Trials')
else
    axis off
end


end