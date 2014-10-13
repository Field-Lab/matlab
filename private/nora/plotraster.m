function raster_plot(xval, fittedGLM)

dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];

figure;
title('Rasters')
rec_rast = xval.rasters.recorded(:,1:bins);
sim_rast = xval.rasters.glm_sim(:,1:bins);

trials   = size(rec_rast,1);
ylim([1 trials*2])

hold on
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    plot(rec1, i_trial, 'k.')
    sim1 = time(find(sim_rast(i_trial,:)));
    if length(sim1) < 4*length(rec1)
        plot(sim1, i_trial + trials, 'r.')
    end
end

xlabel('Time (seconds)')

set(gca,'YTickLabel',{'1','Data',num2str(trials),'Model'},'YTick',[1 trials/2 trials 3*trials/2])

ylabel('Trials')


end