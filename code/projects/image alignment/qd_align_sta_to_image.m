
switch 1
    case 1 % 2007-09-18-4
        % execute image_catalog before this script

        %datarun.piece.photographic_mapping.images.pm2 = imread([server_data_path '2007-09-18-4/slidebook/Image 3']);
        %datarun.piece.photographic_mapping.images.pm32 = imread([server_data_path '2007-09-18-4/slidebook/Image 2']);
        %datarun.piece.photographic_mapping.images.array_edges = imread([server_data_path '2007-09-18-4/slidebook/Image 1']);
        
        scale = 4;

        T_sta_to_monitor = coordinate_transform(datarun,'monitor','scale',scale);

        im = struct;

        im(1).im = f5;%imread([server_data_path '2007-09-18-4/confocal/2009-07-07/2007-09-18-4_63x_stitch_02_relevant/2007-09-18-4_63x_stitch_02_relevant_max.tif']);
        im(1).T_to_array = maketform('composite',TA1,TA1F5);

        im(2).im = a2;
        im(2).T_to_array = TA2;

        %cell_id = 3559; axon_id = 17;  %f3
        %cell_id = 4446; axon_id = 44;  %f4
        %cell_id = 2296; axon_id = 9;  %f2
        %cell_id = 1326; axon_id = 41;  %f4
        %cell_id = 5134; axon_id = 65;  %f4, f5
        cell_id = 413; axon_id = 80;  %f5
        sta_center = rf_com(get_sta(datarun,cell_id));
        im_rad = 5;
        
end



switch 1
    case 1
        % get bounds of RF in monitor coordinates
        rf.im = get_rf(datarun,cell_id,'polarity',true);
        temp = tformfwd(T_sta_to_monitor,[1 1; size(rf.im,2) size(rf.im,1)]);
        rf.im = matrix_scaled_up(rf.im,scale);
        rf.x = temp(:,1);
        rf.y = temp(:,2);
    case 2
        % display monitor alignment image
        rf = load_alignment_image('pm','pm32');
        rf.x = rf.xdata;
        rf.y = rf.ydata;

end

% create transformation from image to monitor
for ii = 1:length(im)
    im(ii).T_image_to_monitor = maketform('composite',...
        fliptform(datarun.piece.T_monitor_to_array),... array to monitor
        im(ii).T_to_array);                           % image to array
end

% identify interesting region of monitor coordinates
temp = tformfwd(T_sta_to_monitor,[sta_center - im_rad; sta_center + im_rad] );
mon_x = temp(:,1)';
mon_y = temp(:,2)';

% transform image
for ii = 1:length(im)
    [im(ii).im_mon,im(ii).mon_x,im(ii).mon_y] = ...
        imtransform(im(ii).im,im(ii).T_image_to_monitor,'xdata',mon_x,'ydata',mon_y,'xyscale',.1);
end

% prepare axon
axon = tforminv(datarun.piece.T_monitor_to_array,axons{axon_id});
axon = traced_cell_points(axon(1:2,:),axon(2:end,:));

% prepare figure
figure(52);clf;num_x = length(im)+1;

% plot image
subplot(1,num_x,1);
imagesc(norm_image(rf.im),'xdata',rf.x,'ydata',rf.y);
hold on
plot(axon(:,1),axon(:,2),'r')
axis image
xlim(mon_x);ylim(mon_y);

% plot each image
for ii = 1:length(im)
    subplot(1,num_x,ii+1);
    imagesc(im(ii).im_mon,'xdata',im(ii).mon_x,'ydata',im(ii).mon_y);
    hold on
    plot(axon(:,1),axon(:,2),'r')
    axis image
end

xlim(mon_x);ylim(mon_y);

% link the axes
linkaxes(get(gcf,'child'))
