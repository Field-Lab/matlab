function plotscatter(xval, fittedGLM)

dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
time     = dt*[1:bins];
rec_rast = xval.rasters.recorded(:,1:bins);
sim_rast = xval.rasters.glm_sim(:,1:bins);
trials   = size(rec_rast,1);

runawaytrial=ones(trials,1);
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    sim1 = time(find(sim_rast(i_trial,:)));
    if length(sim1) > 4*length(rec1)
        runawaytrial(i_trial)=0;
    end
end

runtrials=find(runawaytrial);
convolve=150;
PSTH_rec=zeros(length(runtrials),bins);
PSTH_sim=zeros(length(runtrials),bins);

for i=1:length(runtrials)
    i_trial=runtrials(i);
    PSTH_rec(i,:)=conv(rec_rast(i_trial,:),ones(convolve,1),'same');
    PSTH_sim(i,:)=conv(sim_rast(i_trial,:),ones(convolve,1),'same');
end
MAX=max(max(mean(PSTH_rec)), max(mean(PSTH_sim)))
figure;
set(gcf, 'Position', [100 100 300 700])
subplot(2,1,1)
hold on
c = linspace(1,10,length(mean(PSTH_rec)));
scatter(mean(PSTH_rec),mean(PSTH_sim),5,c);
plot([0 MAX],[0 MAX],'k');
xlim([0 MAX])
ylim([0 MAX])
axis square
xlabel('Recorded Spike Rate')
ylabel('Simulated Spike Rate')
set(gca,'XTick',[],'YTick',[])
subplot(2,1,2)
hold on
scatter(time,mean(PSTH_rec),2,c);
plot(time,mean(PSTH_sim),'k');
xlim([time(1) time(end)])
legend('Recorded','Simulated')

end