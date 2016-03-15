% psth_rather function from code folder with the plotting component removed
function psth_r=psth_raster_noPlotting(mmin,mmax,n,tr,Color,y_scale, counter)


% for i=1:length(tr),
        h=n-tr(1);  
        hh=find(h>=mmin & h<=mmax);
        psth_r=[(h(hh)*1000)', repmat(counter-1,[length(hh),1]);];
% end

hold off;