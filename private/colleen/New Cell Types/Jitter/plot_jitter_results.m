number =[481 813 2103 2163 3288 3816 3888 4280 4937 6336 6512 1534];
sta = zeros(320,640);
figure;
for i = 1:length(number)
    hold on 
    load(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026/Cell ', num2str(number(i)), '.mat'])
    temp = permute(temp, [2,1,3,4]);
    [sig_stixels, params, rf_strength_out] = significant_stixels(temp, 'select', 'thresh', 'thresh', 4.25);
%     plot_sta_(temp)
%     title({'2016-02-17-6 data026' ;['Cell ', num2str(number)]})
%     axis off
%     sta((sig_stixels) == 1) = i;
    sta = sta+ sig_stixels;
%     plot_axes = set_up_fig_or_axes(params.figure);
% image(sta)
% hold on
end

% choose first frame to show
% [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));

% normalize STA color
% sta = norm_image(sta);

% create slider control
% ha = make_loop_slider_list(start_index, 1, size(sta, 4), {@slider_plot, sta, params});

% plot once before any clicks
% slider_plot(ha, [], sta, params);


% function slider_plot(handle, event, sta, params) %#ok<INUSL>
% % display one frame of the STA
% 
% % get the slider position
% ss = round(get(handle,'Value'));
% cla;
% 
% % plot spatial sensitivity
figure; imagesc(sta(10:end-10, 10:end-10));
figure; imshow(sta(10:end-10, 10:end-10));

% image(sta(:,:,:,start_index));
axis image
%set(gca,'xtick',[],'ytick',[])

% if ~isempty(params.overlay)
%     hold on
%     plot(params.overlay(:,1),params.overlay(:,2),'k')
% end

% % title
% title(sprintf('%sframe %d of %d (%d)',params.prefix,ss,size(sta,4),ss-size(sta,4)))





