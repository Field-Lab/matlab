function h3=plot_record_prediction_pc2015_03_09_2(spkCondColl,spkCondCollGLM)

nmov = length(spkCondColl);



xx=cell(nmov,1);
yy=cell(nmov,1);
for imov=1:nmov

%plotraster(x,fittedGLM,'labels',true,'raster_length',24)
[xx{imov},yy{imov}]=plotSpikeRaster(spkCondCollGLM{imov}==1,'PlotType','vertline');
xx{imov}=xx{imov}/(1200);

end



h3= figure;
subplot(2,1,1);
jump=0;
for imov=nmov:-1:1
    
    % Simulated
   
    plot(xx{imov},yy{imov}+jump);
    hold on;
    jump=jump+max(yy{imov});
    
    % Recorded    
    plot(spkCondColl(imov).xPoints/20000,spkCondColl(imov).yPoints+jump)
    hold on
    jump = jump + max(spkCondColl(imov).yPoints);



end
ylim([0,jump]);
xlim([0,max(xx{imov})]);

end