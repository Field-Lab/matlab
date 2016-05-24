%% Figure for EJ's grant on RPE piece stability
% Piece 2016-02-17-1/data000 compared to 2016-02-17-1/data028
% including these two runs, they are 12.27 hours apart
% RGB-8-2-0.48-11111-119.5.xml
% data000: ON parasol cells 1456 1636 1685 1818 2629 2851 3890
% data028: ON parasol cells 1577 1638 1681 1816 2611 2732 2866
clear
close all

dataparam.date = '2016-02-17-1'
dataparam.concatname{1} = 'data000';
dataparam.cell_specification{1} = [1456 1636 1685 1818 2629 2851 3890]; % clear this variable

dataparam.concatname{2} = 'data028';
dataparam.cell_specification{2} = [1577 1638 1681 1816 2611 2732 2866]; % clear this variable

padding = 10;
for num_datarun = 1:size(dataparam.concatname,2);
    clear datarun
    dataparam.file_name = [dataparam.date, '/', dataparam.concatname{num_datarun},'/', dataparam.concatname{num_datarun}];
    
    %% END OF INPUT
    
    datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name, '.neurons'];
    datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name, '.params'];
    datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name, '.sta'];
    
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_all',true);
    opt.load_sta_params.save_rf = 1;
    datarun=load_data(datarun,opt);
    
    cell_indices = get_cell_indices(datarun, dataparam.cell_specification{num_datarun});
    fig_tc{num_datarun} = figure;
    colors = [1 0 0 ; 0 1 0; 0 0 1];
    figure;
        plot_axes = plot_rf_summaries(datarun, dataparam.cell_specification{num_datarun})
title(dataparam.concatname{num_datarun});
        
    for i = 1:length(cell_indices)
        sta = datarun.stas.stas{cell_indices(i)};
        sta = double(sta);
        [~,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
        
        tc(:,1) = datarun.vision.timecourses(cell_indices(i)).r;
        tc(:,2) = datarun.vision.timecourses(cell_indices(i)).g;
        tc(:,3) = datarun.vision.timecourses(cell_indices(i)).b;
        
%         scale = 2;
%         figure;
%         image = norm_image(sta);
%         
%         for j = 1:size(image,3)
%             upsamp_image(:,:,j) =  imresize(image(:,:,j,start_index), scale, 'nearest');
%         end
%         
%         sta_img = imagesc(upsamp_image);
%         
%         hold on
%         
%         the_fit = datarun.stas.fits{cell_indices(i)};
%         ctr = the_fit.mean;
%         rad = the_fit.sd;
%         angle = the_fit.angle;
%         
%         
%         
%         [X,Y] = drawEllipse([ctr*scale rad*scale angle]);
%         hold on
%         plot(X,Y,'Color','k', 'linewidth',2 );
%         axis image
%         hold on
%         set(gca, 'ydir', 'reverse')
%         
%         bounds = autozoom_to_fit(datarun, cell_indices(i), padding, [1,1]);
%         
%         
%         set(gca, 'XLim', bounds(1:2)*scale)
%         set(gca, 'YLim', bounds(3:4)*scale)
%         
%         axis square
%         axis off
%         pix_field_width = datarun.stimulus.field_width*datarun.stimulus.stixel_width;
%         scale_bar_width = 20; %pixels
%         stix_scale_bar = (scale_bar_width/pix_field_width)*datarun.stimulus.field_width;
%         plot([bounds(1)+1; bounds(1)+stix_scale_bar+1]*scale, [bounds(4)-1; bounds(4)-1]*scale, '-k', 'LineWidth', 2)
%         
%         
        
        
        figure(fig_tc{num_datarun});
        hold on
        for c = 1:3
            plot(1/120*2*[0:size(tc,1)-1]', tc(:,c), 'Color', colors(c,:))
        end
        
        title(dataparam.concatname{num_datarun})
    end
    
    
end

