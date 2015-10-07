function h3=plot_record_pred12(spkCondColl,spkCondCollGLM,pred2)
col='rkrkrkrkrkrkrkrkrkr';
nmov = 2;



xx=cell(nmov,1);
yy=cell(nmov,1);

imov=1;
%plotraster(x,fittedGLM,'labels',true,'raster_length',24)
[xx{imov},yy{imov}]=plotSpikeRaster(spkCondCollGLM{1}>0,'PlotType','vertline');
xx{imov}=xx{imov}/(1200);

imov=2;
[xx{imov},yy{imov}]=plotSpikeRaster(pred2{1}>0,'PlotType','vertline');
xx{imov}=xx{imov}/(1200);

h3= figure('Color','w');

jump=0;
icnt=0;
    
    % Simulated
    imov=1;
   icnt=icnt+1;
    plot(xx{imov},yy{imov}+jump,col(icnt));
    hold on;
    jump=jump+max(yy{imov});
    
    % Recorded 
    icnt=icnt+1;
    plot(spkCondColl(imov).xPoints/20000,spkCondColl(imov).yPoints+jump,col(icnt));
    hold on
    jump = jump + max(spkCondColl(imov).yPoints);


    % Simulated
   icnt=icnt+1;imov=2;
    plot(xx{imov},yy{imov}+jump,col(icnt));
    hold on;
    jump=jump+max(yy{imov});

ylim([0,jump]);
xlim([0,max(xx{imov})]);
set(gca,'YTick',[]);
end