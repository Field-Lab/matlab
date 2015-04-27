% AKHeitman
% 2015-01-26

clear; close all;
NSEMmovie = '/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/inputstats_DimFlash_092413Fc12_shift0.mat';
eval(sprintf('load %s', NSEMmovie));
savedir = '/Users/akheitman/NSEM_Home/PrototypePlots/plot_inputNL';



%% Individual Plots of Hinge Mean
inputmean = inputstats.mu_avgIperpix / 256;
xvals_1 = linspace(0,inputmean,100);
xvals_2 = linspace(inputmean,1,100);
ratio_values = .1:.1:1
LW = 4;
MS = 24;
for i_ratio = 1:length(ratio_values)
    rval = ratio_values(i_ratio);
    
    rescale =  1 / (inputmean + rval*(1-inputmean) );
    yvals_1 = rescale*xvals_1;
    yvals_2 = rescale * rval * (xvals_2 - inputmean) + rescale *inputmean;
    
    clf; hold on;
    plot([xvals_1,xvals_2],[yvals_1,yvals_2], 'linewidth',LW); 
    plot( inputmean,rescale*inputmean, 'k.', 'markersize', MS);
    xlim([0,1]); ylim([0,1]);
    %text(inputmean,0,'Mean of Input');
    %text(inputmean,rescale*inputmean,' \leftarrow Hinge At Input Mean','FontSize',18)
    
    title(sprintf('HingeAboutMean_IncrementDecrementRatio_%1.1e', rval), 'interpreter','none');
    xlabel('Pre-NonLinearity');
    ylabel('Post-NonLinearity');
    
    
    orient landscape
    eval(sprintf('print -dpdf %s/HingeMean_Inc_%1.1e.pdf', savedir, rval)); 
end

%%
offpar = [1.5000 1 1 0.5000];

onpar_2824 = [2.5000 1 0 0.5000];
onpar_6799 = [2.5000 1 0.5000 0];
onpar_3167 = [2.5000 1 0 0.5000];    
onpar_5660 = [1.5000 1 0.5000 1];   


oldstim = linspace(0,1,1000);
quartile_1 =   1:250;
quartile_2 = 251:500;
quartile_3 = 501:750;
quartile_4 = 751:1000;


LW = 4;
MS = 24;
for i_plot = 1:5
    newstim = oldstim;
    if i_plot == 1, slope = offpar; cellid = 'offpar'; end
    if i_plot == 2, slope = onpar_2824; cellid = 'onpar_2824';end
    if i_plot == 3, slope = onpar_6799; cellid = 'onpar_6799';end
    if i_plot == 4, slope = onpar_3167;cellid = 'onpar_3167'; end
    if i_plot == 5, slope = onpar_5660; cellid = 'onpar_5660'; end
        
    offset1 = 0;
    offset2 = .25*slope(1);
    offset3 = .25*(slope(1)+slope(2));
    offset4 = .25*(slope(1)+slope(2)+slope(3));

    newstim(quartile_1) = slope(1) * (newstim(quartile_1) - .00) + offset1;
    newstim(quartile_2) = slope(2) * (newstim(quartile_2) - .25) + offset2;
    newstim(quartile_3) = slope(3) * (newstim(quartile_3) - .50) + offset3;
    newstim(quartile_4) = slope(4) * (newstim(quartile_4) - .75) + offset4;

    clf;plot(linspace(0,1,1000), newstim,'linewidth',LW);hold on;
    xlabel('Pre-NonLinearity');
    ylabel('Post-NonLinearity');
    
    plot(oldstim(250),newstim(250),'k.','markersize', MS);
    plot(oldstim(500),newstim(500),'k.','markersize', MS);
    plot(oldstim(750),newstim(750),'k.','markersize', MS);
    
    %text(oldstim(250),newstim(250),' \leftarrow Hinge At Input Val .25','FontSize',18);
    %text(oldstim(500),newstim(500),' \leftarrow Hinge At Input Val .5','FontSize',18);
    %text(oldstim(750),newstim(750),' \leftarrow Hinge At Input Val .75','FontSize',18);
    
    title(sprintf('PieceLin_foursegment_%s: Hinge at .25, .5, .75', cellid), 'interpreter','none');
    orient landscape
    eval(sprintf('print -dpdf %s/PieceLin_foursegment_%s.pdf', savedir, cellid)); 
end



