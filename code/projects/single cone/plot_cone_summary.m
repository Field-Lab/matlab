function plot_cone_summary(datarun, varargin)
% plot_cone_summary     Plot cone mosaic, ROI, and DRP
%
% usage:  plot_cone_summary(datarun, varargin)
%
% arguments:     datarun - datarun struct
%               varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% num_bins          20     	how many bins for the DRP
% bin_width        	1       bin width
% fig_drp          	0      	which figure to plot mosaic, ROI, and DRP.  if 0, make new figure.
% fig_class        	0      	which figure to plot classification
% cone_roi         	[]      binary vector of cones in the ROI.  If empty, not plotted.
%
%
%
% 2008-10 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('num_bins', 20);
p.addParamValue('bin_width', 1);
p.addParamValue('fig_drp', 0);
p.addParamValue('fig_class', 0);
p.addParamValue('cone_roi', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% SET UP FIGURES


% clear figure (or make new)
set_up_fig_or_axes(params.fig_drp);

% make correct units
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize=28;
xLeft = (22-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])

% make plot axes
drp_plot_axes = subplot_axes(params.fig_drp,[.07 .05 .93 .9],.25,.15,2,3);



% clear figure (or make new)
set_up_fig_or_axes(params.fig_class);

% make correct units
set(gcf,'PaperUnits','centimeters')
xSize = 20; ySize=28;
xLeft = (22-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])

% make plot axes
class_plot_axes = subplot_axes(params.fig_class,[.07 .05 .93 .9],.25,.15,2,2);





% put variables into nicer format
centers = datarun.cones.centers;
Lc = centers(datarun.cones.types=='L',:);
Mc = centers(datarun.cones.types=='M',:);
Sc = centers(datarun.cones.types=='S',:);
Uc = centers(datarun.cones.types=='U',:);




%%%%%%%%%%%%%%%%%%%%%  cone mosaic  %%%%%%%%%%%%%%%%%%%%%

plot_cone_mosaic(datarun,'fig_or_axes',drp_plot_axes{1},'cone_size',2.8)
set(drp_plot_axes{1},'XTick',[],'YTick',[])
title(drp_plot_axes{1},printable_name(datarun))





%%%%%%%%%%%%%%%%%%%%%  DRP  %%%%%%%%%%%%%%%%%%%%%

% parameters
num_bins = params.num_bins;
bin_width = params.bin_width;

[drp,bin_centers,extras] = density_recovery_profile(centers,num_bins,bin_width);
axes(drp_plot_axes{3})
bar(bin_centers,drp)
xlabel('distance (pixels)');ylabel('density (cones/pixel^2)')

% add mean density, and effective radius
hold on
xlim = get(gca,'XLim');
plot(xlim,[1 1]*extras.density,'r',[1 1]*extras.eff_rad,[0 extras.density],'r')
text(xlim(2)*0.9,extras.density*1.1,sprintf('mean density: %0.3f\nradius: %0.1f',extras.density,extras.eff_rad),...
    'HorizontalAlignment','Right','VerticalAlignment','Bottom')

% set ylim
best_ylim = [0 max(drp)*1.3];
set(gca,'YLim',best_ylim)




%%%%%%%%%%%%%%%%%%%%%  cone ROI  %%%%%%%%%%%%%%%%%%%%%

if ~isempty(params.cone_roi)
    % plot it
    plot_cone_mosaic(datarun,'fig_or_axes',drp_plot_axes{2},'cone_size',2.8,'cone_roi',params.cone_roi)
    set(drp_plot_axes{2},'XTick',[],'YTick',[])
else
    % clear if unused
    cla(drp_plot_axes{2})
end





%%%%%%%%%%%%%%%%%%%%%  DRP, grouped by cone type  %%%%%%%%%%%%%%%%%%%%%

% from L cones
[drp_L,bin_centers] = density_recovery_profile(Lc,num_bins,bin_width);
drp_M = density_recovery_profile(Mc,num_bins,bin_width,'reference_centers',Lc);
drp_S = density_recovery_profile(Sc,num_bins,bin_width,'reference_centers',Lc);

data = [drp_S; drp_M; drp_L]';

axes(drp_plot_axes{5})
bar(bin_centers(1:end),data(1:end,:),'stacked')
set(gca,'YLim',best_ylim)
xlabel('distance from L cones (pixels)');ylabel('density (cones/pixel^2)')
d1=data;


% from M cones
[drp_L,bin_centers] = density_recovery_profile(Lc,num_bins,bin_width,'reference_centers',Mc);
drp_M = density_recovery_profile(Mc,num_bins,bin_width);
drp_S = density_recovery_profile(Sc,num_bins,bin_width,'reference_centers',Mc);

data = [drp_S; drp_M; drp_L]';

axes(drp_plot_axes{6})
bar(bin_centers(1:end),data(1:end,:),'stacked')
set(gca,'YLim',best_ylim)
xlabel('distance from M cones (pixels)');ylabel('density (cones/pixel^2)')
d2=data;



% from S cones
[drp_L,bin_centers] = density_recovery_profile(Lc,num_bins,bin_width,'reference_centers',Sc);
drp_M = density_recovery_profile(Mc,num_bins,bin_width,'reference_centers',Sc);
drp_S = density_recovery_profile(Sc,num_bins,bin_width);

data = [drp_S; drp_M; drp_L]';

axes(drp_plot_axes{4})
bar(bin_centers(1:end),data(1:end,:),'stacked')
set(gca,'YLim',best_ylim)
xlabel('distance from S cones (pixels)');ylabel('density (cones/pixel^2)')
d3=data;



%dsum=0.475*d1+0.475*d2+0.05*d3;
%dsum=d1+d2+d3;
%axes(drp_plot_axes{2})
%bar(bin_centers(1:end),dsum(1:end,:),'stacked')







%%%%%%%%%%%%%%%%%%%%%  cone classification  %%%%%%%%%%%%%%%%%%%%%


% plot cone classification
plot_cone_classification(datarun,'foa_3d',class_plot_axes{1},'foa_2d',class_plot_axes{2},...
    'foa_pie',class_plot_axes{3},'foa_hist',class_plot_axes{4},'dot_size',2);

% add title
title(class_plot_axes{1},printable_name(datarun))

% clear unused panel
%set(class_plot_axes{4},'visible','off')


