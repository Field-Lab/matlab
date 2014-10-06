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
colorbar