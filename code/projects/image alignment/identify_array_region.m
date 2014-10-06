% identify in an original image where the array is located
% assumes a transformation from image to array has already been made;

% parameters
array_type = 512;
orig_image = f1;
tform = maketform('composite',TA1,TA1F1);


% get positions of transformed electrode corners
ep = electrode_positions(array_type);
switch array_type
    case 512
        corner_electrodes = [512 129 249 392 512];
end
corner_positions = ep(corner_electrodes,:);
% apply transformation
corner_positions = tforminv(tform,corner_positions);


% plot original image
figure;imagesc(orig_image);colormap gray;axis image;hold on;
% plot transformed lcoations
plot(corner_positions(:,1),corner_positions(:,2),'r')
% save to desktop
print(gcf,'/snle/home/gauthier2/Desktop/temp','-dpdf')
