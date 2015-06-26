function slider_plot(handle, event, sta) %#ok<INUSL>
% display one frame of the STA

% get the slider position
ss = round(get(handle,'Value'));
cla;
% plot spatial sensitivity


imagesc(sta(:,:,ss));
colormap gray
caxis([min(sta(:)),max(sta(:))]);
colorbar

%set(gca,'xtick',[],'ytick',[])


% title
title(sprintf('frame %d of %d (%d)',ss,size(sta,3),ss-size(sta,3)))
