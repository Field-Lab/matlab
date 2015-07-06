function plot_axes = plot_zoomed_STAs(datarun, cell_specification, params)    
params.foa = 0;
params.clear= 0;
plot_axes = set_up_fig_or_axes(params.foa,params.clear);

[cell_indices, cell_type_name] = get_cell_indices(datarun, cell_specification);
padding = params.padding;
for i = 1:1%length(cell_indices)
    axis equal
        bounds = autozoom_to_fit(datarun, cell_indices(i), padding, [1,1]);
    
    
    sta=squeeze(datarun.stas.stas{cell_indices(i)});
    if size(sta, 3) ~= 3
        [~,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));
        image('Parent',plot_axes,'CData',sta(:,:,start_index))
        
    else
        [~,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
        
        sta = norm_image(sta);
        image('Parent',plot_axes,'CData',sta(:,:,:,start_index))
    end
    
%     if vision ~= 1
%         the_fit.mean  = [datarun2.matlab.sta_fits{cell_numbers(k)}.center_point_x  datarun2.matlab.sta_fits{cell_numbers(k)}.center_point_y];
%         the_fit.sd  = [datarun2.matlab.sta_fits{cell_numbers(k)}.center_sd_x  datarun2.matlab.sta_fits{cell_numbers(k)}.center_sd_y];
%         the_fit.angle = datarun2.matlab.sta_fits{cell_numbers(k)}.center_rotation_angle;
%     else
        the_fit = datarun.stas.fits{cell_indices(i)};
%     end
    
    % get points of an ellipse with these parameters
    [X,Y] = drawEllipse([the_fit.mean the_fit.sd the_fit.angle]);
    hold on
    plot( X,Y,'Color','k', 'LineWidth',1, 'Parent',plot_axes)

    
    
    set(plot_axes, 'XLim', bounds(1:2))
        set(plot_axes, 'YLim', bounds(3:4))

    axis square
    axis off
    pix_field_width = datarun.stimulus.field_width*datarun.stimulus.stixel_width;
    scale_bar_width = 15; %pixels
    stix_scale_bar = (scale_bar_width/pix_field_width)*datarun.stimulus.field_width;
    plot([bounds(1)+1; bounds(1)+stix_scale_bar+1], [bounds(4)-1; bounds(4)-1], '-k', 'LineWidth', 2, 'Parent', plot_axes)
       
   
end
