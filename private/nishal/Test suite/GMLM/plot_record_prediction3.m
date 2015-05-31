function h3=plot_record_prediction3(spkCondColl,spkCondCollGLM,pred2)
col='rkmrkmrkmrkmrkmrkmrkmrkmrkmrk';
nmov = length(spkCondColl);



xx=cell(nmov,1);
yy=cell(nmov,1);
for imov=1:nmov

%plotraster(x,fittedGLM,'labels',true,'raster_length',24)
[xx{imov},yy{imov}]=plotSpikeRaster(spkCondCollGLM{imov}>0,'PlotType','vertline');
xx{imov}=xx{imov}/(1200);

end

for imov=1:nmov

%plotraster(x,fittedGLM,'labels',true,'raster_length',24)
[xx2{imov},yy2{imov}]=plotSpikeRaster(pred2{1}>0,'PlotType','vertline');
xx2{imov}=xx2{imov}/(1200);

end



h3= figure;

jump=0;
icnt=0;
for imov=nmov:-1:1
     % Simulated2
   icnt=icnt+1;
    plot(xx2{imov},yy2{imov}+jump,col(icnt));
    hold on;
    jump=jump+max(yy2{imov}); 
    
    % Simulated
   icnt=icnt+1;
    plot(xx{imov},yy{imov}+jump,col(icnt));
    hold on;
    jump=jump+max(yy{imov});
    
    % Recorded 1
    icnt=icnt+1;
    plot(spkCondColl(imov).xPoints/20000,spkCondColl(imov).yPoints+jump,col(icnt));
    hold on
    jump = jump + max(spkCondColl(imov).yPoints);

end
ylim([0,jump]);
xlim([0,max(xx{imov})]);

end