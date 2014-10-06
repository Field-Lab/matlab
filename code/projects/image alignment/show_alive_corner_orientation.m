function show_alive_corner_orientation(array_id, rig)
% SHOW_ALIVE_CORNER_ORIENTATION    Help for aligning alive image to array
% 
% The fast method for aligning live image to array is to click the corner
% points.  The question is what order to click the points in.  This answers
% that question.
%
% 2010-05 phli
%

old_hold = ishold;


array_info = load_array_info(array_id);

pos = array_info.positions;
im_pos = tformfwd(array_info.T_array_to_array_image, pos);

clf;
hold on;
set(gca, 'Visible', 'off');
set(gcf, 'Color', 'w');
switch rig
    case 'C'
        axis xy; % Luckily this appears to be an okay way to handle this...
    otherwise
        axis ij;
end
plot(im_pos(:,1), im_pos(:,2), 'go');
axis equal;

corner_electrodes = array_info.corner_electrodes;
plot(im_pos(corner_electrodes,1), im_pos(corner_electrodes,2), 'bs');
for i = 1:length(corner_electrodes)
    elec = corner_electrodes(i);
    text(im_pos(elec,1), im_pos(elec,2), ['   ' num2str(i) ' (e ' num2str(elec) ')'], 'Color', 'b');
end


if ~old_hold
    hold off;
end