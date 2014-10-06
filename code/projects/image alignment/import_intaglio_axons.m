function [axons T] = import_intaglio_axons(svg_file_path,images,T_align)
% import_intaglio_axons     import axons paths from an SVG file, and transform to array coordinates
%
%
%      SEE GRAPHICAL DOCUMENTATION IN FILE "import_intaglio_axons.intaglio"
%
%
% usage:  import_intaglio_axons(svg_file_path,images,T_align)
%
% arguments:     svg_file_path - path to SVG file
%                                   the file should be exported from intaglio
%                                   the "images" layer should contain images which have known landmark locations
%                                       usually the axons are actually visible in these images
%                                   the "axons" layer should contains paths which trace out the axon trajectories
%
%                       images - cell array describing the images in the SVG file
%                                   images{jj}    - information about image jj, an "axon image"
%                                   images{jj}{1} - name of the image in the SVG file, e.g. '18f20030.tiff'
%                                             {2} - name of the image in matlab, e.g. 'f2'
%                                             {3} - transformation from array coordinates to image coordinates 
%                                             {4} - coordinates of landmarks in the alignment image
%                                             {5} - coordinates of landmarks in the axon image
%
%                      T_align - transformation from the alignment image to array coordinates
%
%
% outputs:     axons - cell array of axon paths, each an Nx2 matrix
%                       row 1 is the center of the soma
%                       row 2 is where the axon leaves the soma
%                       subsequent rows give the rest of the axon path
%                       the order of axons in the cell array is based on sorting the x-coordinate of the soma centers
%
%
%
% 2010-03 gauthier
%






% LOAD SVG FILE


% load file
xDoc = xmlread(svg_file_path);

% find which layers contain the axons and images

% look through each layer
for ii = 0:xDoc.getElementsByTagName('svg').item(0).getLength - 1
    % if its name matches any expected name, note its number
    switch xDoc.getElementsByTagName('svg').item(0).item(ii).getAttribute('id').toCharArray'
        case 'axons'
            axon_layer = ii;
        case 'images'
            image_layer = ii;
    end
end
if ~exist('axon_layer','var'); error('could not identify which layer contains axons');end
if ~exist('image_layer','var'); error('could not identify which layer contains images');end




% GET AXONS IN SVG COORDINATES


% initialize axon storage variable
axons_svg = cell(0);

% note number of items in this layer (some are axons, some are not)
num_items = xDoc.getElementsByTagName('svg').item(0).item(axon_layer).getLength;

% loop through each item
for aa = 0:num_items-1

    % get the item
    this_item = xDoc.getElementsByTagName('svg').item(0).item(axon_layer).item(aa);

    % identify whether it is an axon (i.e. a path object with no fill)
    if strcmp(this_item.getNodeName,'path') && strcmp(this_item.getAttribute('fill').toCharArray','none')

        % get text of axon coordinates
        char_points=this_item.getAttribute('d').toCharArray;

        % convert to numbers
        axons_svg{length(axons_svg)+1} = svg_path_to_coordinates(char_points);

    end
end

% plot
if 0
    figure;hold on
    for aa=1:length(axons_svg)
        plot(axons_svg{aa}(:,1),axons_svg{aa}(:,2),'Color',rand(1,3))
        pause
    end
end




% COMPUTE TRANSFORM FROM SVG COORDINATES TO ELECTRODE COORDINATES


% identify which objects in the "images" layer are actually images

% initialize list of images
image_indices = [];

% go through each object...
for ii = 0:xDoc.getElementsByTagName('svg').item(0).item(image_layer).getLength-1
    % if it is an image
    if ~isempty(regexp(xDoc.getElementsByTagName('svg').item(0).item(image_layer).item(ii).toString.toCharArray','image','once'))
        %note its index
        image_indices = [image_indices ii];
    end
end

% note the total number of images
num_images = length(image_indices);


% initialize storage of alignment points from each image
all_alignment_points = [];
all_svg_points = [];


% go through each image
for ii = 1:num_images

    % get the image
    this_image = xDoc.getElementsByTagName('svg').item(0).item(image_layer).item(image_indices(ii));
    
    % find which axon image this is
    for jj=1:length(images)
        if strcmp(images{jj}{1},this_image.getAttribute('xlink:href').toCharArray')
            % get the cell body locations that were used to align this fixed image with the alive image
            axon_points = images{jj}{4};
            alignment_points = images{jj}{5};
            break
        end
    end
    
    % check to be sure it was found
    if ~exist('alignment_points','var')
        error('the svg file image ''%s'' was not found in provided list of images.',this_image.getAttribute('xlink:href').toCharArray')
    end
    
    % make transform from image coordinates to SVG coordinates 
    
    % http://www.w3.org/TR/SVG11/coords.html#TransformMatrixDefined
    
    % load matrix transform text
    matrix_transform = this_image.getAttribute('transform').toCharArray';
    % extract numbers
    matrix_transform = textscan(matrix_transform(8:end-1),'%f','delimiter',',');
    % create matlab transform, from image coordinates to SVG coordinates
    T_axon_to_svg = maketform('projective',[reshape(matrix_transform{1},2,3); 0 0 1]');

    
    % add alive image points to list
    all_alignment_points = [all_alignment_points; alignment_points];
    
    % add fixed image points to list, after transforming them to SVG file coordinates
    all_svg_points = [all_svg_points; tformfwd(T_axon_to_svg,axon_points)];
    
    % clean up 
    clear alignment_points axon_points

end

% transform alignment points to electrode coordinates
all_elec_points = tformfwd(T_align,all_alignment_points);

% compute transformation to electrode coordinates
T = cp2tform(all_elec_points,all_svg_points,'lwm');



% APPLY TRANSFORMATION TO AXONS
    
% for each axon
for aa=1:length(axons_svg)
    % transform to array coordinates
    ax_array{aa} = tforminv(T,axons_svg{aa});
end



% SORT THE ORDER OF THE AXONS

% initialize axon starting point storage
axon_starts = [];

% get starting point of each axon
for aa=1:length(ax_array)
    axon_starts = [axon_starts; ax_array{aa}(1,1)];
end

% identify order 
[junk,new_order] = sort(axon_starts);

% re-order
axons = {ax_array{new_order}};


% show axon count
fprintf('loaded %d axons from %s\n\n',length(axons),svg_file_path)
