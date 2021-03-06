function [uop_bps crm_bps] = rastbps_comp(raster,t_bin,psfilter,options)
% STAND-ALONE
% This program does xyz

%%% PURPOSE %%%
% Compute a Conditioned Rate Model (CRM) from the raster and psfilter
% Compute the Uncoditioned Optimal Rate Model from raster alone
% Compute corresponding Bits_Per_Spike Code
%
%%% NOTES %%%
% Implement 10 msec guassian smoothing time for all spikes
% Uses post-spike filter generated previously (glm fits) as means for
% conditioning
% Fits run through GLM like fit procedure
% Base Model features only 3 free parameters
%
%%% INPUTS  %%%
% raster: rows are repitions, columns times, binary 0 1 for spike 
% t_bin: time duration of each bin in seconds
% psfilter: log(gain change) induced by a spike
%           concept directly from the GLM
%
%%% OUTPUTS %%%
% uop_bps: Bits per Spike of the Unconditioned Optimal Model
% crm_bps: 
%
%%% OUTSIDE CALLS %%%
% NONE
%
%%% KEY MATLAB CALLS %%%
% fminunc
%
% AKHEITMAN 2015-03-01
% Version 0_0  confirmed to work 2015-03-01
% Version 0_1  Full standalone code: confirmated to work 2015-03-01

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% TEST CODE %%%
%{
% INPUT 
load rastbps_comp_test.mat
[uop_bps crm_bps] = rastbps_comp(raster,t_bin,psfilter);
display(sprintf('uop_bps = %d; crm_bps = %d', uop_bps, crm_bps))

if uop_bps == uop_bps && crm_bps = crm
% DESIRED OUTPUT
% uop_bps = 2.895830e-01; crm_bps = 4.697288e-01
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOTTING COMMANDS
&&&&&&&&&&&&&&&&&&
&&&&&&
%%% PLOTTING COMMANDS %%%%%
%s,'interpreter','none')  to make certain symbols act normally!
axis square

figname = sprintf('%s/%s', plotdir,name_figure);
saveas(gcf,figname,'fig')
%eval(sprintf('print -dfig %s/%s.fig',plotdir,name_figure));
eval(sprintf('print -dpdf %s/%s.pdf',plotdir,name_figure));
eval(sprintf('print -depsc %s/%s.eps',plotdir,name_figure));  
% Troubleshooting template
if exist('troubleshoot','var') && troubleshoot.doit
    clf
    subplot(3,2,[3 5]); hist(double(fitmovie(:)),20); set(gca, 'fontsize', 10); title('histogram of full raw stim')
    subplot(3,2,[4 6]); hist(           stim(:) ,20); set(gca, 'fontsize', 10); title('histogram of normed stim over ROI') 
    subplot(3,2,[1 2]);   set(gca, 'fontsize', 10); axis off
    c = 0;
    c=c+1; text(-.1, 1-0.1*c,sprintf('Trouble shooting: %s',troubleshoot.name));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Specifically: Stimulus Normalizing Component of "glm-nospace"' ));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
    c=c+1; text(-.1, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')) );
    
    orient landscape
    eval(sprintf('print -dpdf %s/%s_stimnorm_glmnospace.pdf', troubleshoot.plotdir, troubleshoot.name));
end 

[n1] = hist(values, hist_x);
bar(hist_x,n1,plotstring);
xlim([xmin,xmax]);
 c = 0; 
    c = c+2; text(-.1, 1-0.1*c,sprintf('Metric: %s, Stimulus: %s', raster_computation, stimtype),'interpreter','none');

raster_selfprediction.timestamp           = datestr(clock);
raster_selfprediction.mfile_name          = mfilename('fullpath');


% lots of axis manipulation    
subplot(3,2,3)  % label in 50 msec intervals
set(gca, 'fontsize', 12)
dur         = 1000*frames*tstim;
msec_tick   = 50:50:dur;
frame_tick  = floor(msec_tick / (1000*tstim));
pixel_tick  = 3:3:space;
col_ax = [min(K(:)),max(K(:))];
imagesc(Kxt',col_ax);
set(gca,'xtick',pixel_tick);
xlabel('Space [pixel]')
set(gca, 'ytick', frame_tick, 'yticklabel',msec_tick)
ylabel('Time [msec]')
title('Stimulus Filter 1D-Space by Time')  
colorbar
subplot(3,2,4)
set(gca, 'fontsize', 12)
imagesc(reshape(spfilter,[space,space])); colormap jet
axis image
set(gca,'xtick',pixel_tick);
set(gca,'ytick',pixel_tick);
xlabel('Space [pixel]')
ylabel('Space [pixel]')
mx_time = round(1000*tstim*mx_fr);
title(sprintf('Rk1 Spatial Component'))
colorbar; colorbar('fontsize',12)


%% Plotting Function template
plot_Input_NL_Distribution(changes_cell,plotversion,load_exps,avg_exps)
base_figuredir  = '/Users/akheitman/NSEM_Home/PrototypePlots/plot_inputNL';
figuredir       = sprintf('%s/plot_Input_NL_Distribution_PLOTVERSION_%s',base_figuredir,plotversion);
if ~exist(figuredir, 'dir'), mkdir(figuredir); end
plotNL_v0(figuredir)

function plotNL_v0(data, plotdir)
LW2 = 1.5; LW = 12;
figure(1); clf; hold on;
figname = sprintf('%s/%s', plotdir,name_figure);
saveas(gcf,figname,'fig')
eval(sprintf('print -dpdf %s/%s.pdf',plotdir,name_figure));
eval(sprintf('print -depsc %s/%s.eps',plotdir,name_figure));   
end

