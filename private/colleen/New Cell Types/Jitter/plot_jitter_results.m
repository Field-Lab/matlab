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



[~,git_hash_string] = system('git rev-parse HEAD');
fprintf('Git Hash: %s \n', git_hash_string);


%% Mosiac from Vision
dataparam.date='2016-02-17-6';
dataparam.concatname='data025';
dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
% dataparam.save_path = ['/Volumes/Lab/Users/crhoades/JitterMovie/', dataparam.date, '/', dataparam.concatname, '/'];
select_cells = 1;
if select_cells == 1
dataparam.cell_specification = [482 813 1537 2103 2167 3288 3694 3889 4326 4939 6336 6517];
    
end
dataparam.cell_type = {'all'};

% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

cell_ids = get_cell_indices(datarun, dataparam.cell_specification)

for i = 1:length(cell_ids)
datarun.vision.stas.sta{cell_ids(i)}
end
