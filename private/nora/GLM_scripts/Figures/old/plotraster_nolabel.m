function plotraster_nolabel(xval, fittedGLM)

dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];

rec_rast = xval.rasters.recorded(:,1:bins);
sim_rast = xval.rasters.glm_sim(:,1:bins);

trials   = size(rec_rast,1);

hFig1=figure(1);
set(hFig1, 'Position', [100 100 500 250])
hold on
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'k.')
end

axis off

hFig2=figure(2);
set(hFig2, 'Position', [100 100 500 250])
hold on
for i_trial = 1:trials
    sim1 = time(find(sim_rast(i_trial,:)));
    if length(sim1) < 4*length(rec1)
        plot(sim1, i_trial, 'r.')
    end
end

axis off


end