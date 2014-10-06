% align electrode positions to an image that shows the array
%
% 1) enter parameters
% 2) execute script
% 3) follow instructions printed out in matlab command window
% 4) save the results using the command:
%       save(['/snle/lab/Experiments/Array/Analysis/PIECE/images/' T_name],T_name,'array_corners','alignment_locs')
%
% 2009  gauthier
%



% PARAMETERS

im_name = a1; % which image
T_name = 'TA1'; % name of variable for storing the transformation
array_type = 512; % array type (512, 519, 61)
fig = 1; % which figure to plot in
plot_rad = 200; % radius in pixels for zooming in on electrodes






% choose array type
switch array_type
    case 512
        instruction_text = ['for image in Rig B, start at upper left, go clockwise\n' ...
            'at "blunt" corners, click the electrode which is furthest from the counter-clockwise corner\n'];
        ep = electrode_positions(512);
        electrode_list = [129 249 392 512];
    otherwise
        error('not supported yet!')
end
num_points = length(electrode_list);

% get locations of corner electrodes
alignment_locs = ep(electrode_list,:);


% do it once to get approximate location
figure(fig);clf;imagesc(im_name);axis image;
fprintf('click on corner electrodes.  these clicks only need to be approximately correct\n%s',instruction_text);disp(electrode_list)
rough_points = zeros(num_points,2);
for pp=1:num_points
    waitforbuttonpress
    fum = get(gca,'CurrentPoint');
    rough_points(pp,:) = fum(1,1:2); % button down detected
end

% do it again, this time zoomed in on the electrode
fprintf('click on corner electrodes again (this time the window will be zoomed in on the electrodes)\n')
array_corners = zeros(num_points,2);
for pp=1:num_points
    % plot on zoomed in region
    figure(fig);clf;imagesc(im_name);axis image;
    set(gca,'xlim',plot_rad*[-1 1]+rough_points(pp,1),...
        'ylim',plot_rad*[-1 1]+rough_points(pp,2))
    waitforbuttonpress
    fum = get(gca,'CurrentPoint');
    array_corners(pp,:) = fum(1,1:2); % button down detected
end
close(fig);

% make array transformation
T_temp = cp2tform(array_corners, alignment_locs, 'projective');
eval([T_name '=T_temp;'])



% create transformation from fixed image to alive image
if 0
    cpselect(f1,a1)
    TA1F1 = cp2tform(input_points,base_points,'lwm');
    save('/snle/lab/Experiments/Array/Analysis/PIECE/images/TA1F1','TA1F1','input_points','base_points')
end