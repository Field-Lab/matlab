function [spkCondCollGLM,h]=plot_GLM_prediction_movie_scaled_pc2015_03_09_2(cellID,condMov,GLMPath,mov_scales)

col='krkrkrkrkrkrkr';

nTrials=30;
load(sprintf(strcat(GLMPath,'CellID_%d.mat'),cellID));
nmov =length(condMov);

spkCondCollGLM=cell(nmov,1);
xx=cell(nmov,1);
yy=cell(nmov,1);
for imov=1:nmov
    fittedGLM_scaled = fittedGLM;
    fittedGLM_scaled.linearfilters.Stimulus.Filter=fittedGLM_scaled.linearfilters.Stimulus.Filter*mov_scales(imov);
    fittedGLM_scaled.linearfilters.Stimulus.space_rk1 = fittedGLM_scaled.linearfilters.Stimulus.space_rk1*mov_scales(imov);
x=GLM_predict(fittedGLM_scaled, condMov{imov}, nTrials);
%plotraster(x,fittedGLM,'labels',true,'raster_length',24)
[xx{imov},yy{imov}]=plotSpikeRaster(x.rasters.glm_sim==1,'PlotType','vertline');
xx{imov}=xx{imov}*x.rasters.bintime;

spkCondCollGLM{imov} = x.rasters.glm_sim;
end

h= figure;
subplot(2,1,1);
jump=0;
for imov=nmov:-1:1
plot(xx{imov},yy{imov}+jump,col(imov));
hold on
jump = jump + max(yy{imov});
end
ylim([0,jump]);
xlim([0,max(xx{imov})]);

end