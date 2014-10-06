function [left_image, right_image] = qd_compute_3d_projection(I)

% compute straigt on projection (left_image) and projection rotated by 10
% degrees along the vertical axis (right_image).  both are normalized to
% have a maximum of 1.
% assume image data is a 4D matrix where the third dimension is color (R,G,
% and B) and the fourth dimension is the z position, i.e. depth.


y_slices = size(I,1);  %number of slices in the horizontal dimension
x_slices = size(I,2);  %number of slices in the vertical dimension
z_slices = size(I,4);  %depth slices

theta = 10;

edge_trim_size = ceil((1-cosd(theta))*y_slices);

clear result
clear z_section
clear z_slice
clear temp

tic

%go through each row (x_slice) of the image
for xx = 1:x_slices
    
    % compute the projection, rotated by theta
    for cc = 1:3
        result(xx,:,cc) = radon(reshape(I(xx,:,cc,:),y_slices,z_slices)',theta);
    end

    
    % previous algorithms, not used b/c they're a little slower or less elegant
    
%     % 1) reconstruct the z section into its own matrix
%     for zz = 1:z_slices
%         z_section(zz,:,1:3) = I(xx,:,1:3,zz);
%     end
%     for cc = 1:3
%         result(xx,:,cc) = radon(z_section(:,:,cc),theta);
%     end
% 
%     % 2) more direct
%     for cc = 1:3
%         z_slice(:,:) = I(xx,:,cc,:);
%         result(xx,:,cc) = radon(z_slice',theta);
%     end
    
end

toc


result = result(:,1 + edge_trim_size:y_slices - edge_trim_size,:);


%normalize
right_image = result/max(max(max(result)));

%use the mean to compute the straight on projection
left_image = mean(I,4);
left_image = left_image/max(max(max(left_image)));
