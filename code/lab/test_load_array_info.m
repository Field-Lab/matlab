% this code tests the output of the function load_array_info.
% for each kind of array, plot the array image with electrodes, array abounds, and cornes overlaid.


array_id_list = [100 500 1530 1512];
    

for array_id = array_id_list

    % get array info
    array_info = load_array_info(array_id);
    
    % plot image
    figure(array_id);clf;
    imagesc(array_info.image);axis image;hold on
    
    % get data in array_image coordinates
    pos = tformfwd(array_info.T_array_to_array_image,array_info.positions);
    corners = tformfwd(array_info.T_array_to_array_image,array_info.corners);
    bounds = tformfwd(array_info.T_array_to_array_image,[array_info.x_bounds' array_info.y_bounds']);
    
    % plot electrodes with numbers
    plot(pos(:,1),pos(:,2),'.')
    
    % with numbers
    for ee=1:size(array_info.positions,1)
            text(pos(ee,1),pos(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',12,...
            'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
    
    % plot corners
    plot(corners([1:end 1],1),corners([1:end 1],2),'m')
    
    % plot x, y bounds
    xx = bounds(:,1);
    yy = bounds(:,2);
    plot([xx(1) xx(2) xx(2) xx(1) xx(1)],[yy(1) yy(1) yy(2) yy(2) yy(1)],'k')
    
    % title with shape, spacing
    title(sprintf('array shape: ''%s'', array spacing: %d',array_info.shape,array_info.spacing))
    
end

