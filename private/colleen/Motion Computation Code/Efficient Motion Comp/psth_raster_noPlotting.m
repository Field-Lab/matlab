% psth_rather function from code folder with the plotting component removed
function psth_r=psth_raster_noPlotting(mmin,mmax,n,tr,Color,y_scale)
psth_r=[];

for i=1:length(tr),
        h=n-tr(i);  
        hh=find(h>=mmin & h<=mmax);
        psth_r=[psth_r; (h(hh)*1000)', repmat(length(tr)-i,[length(hh),1]);];
end

hold off;