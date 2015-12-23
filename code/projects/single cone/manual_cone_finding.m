function [new_locations, clicked_stixels] = manual_cone_finding(datarun, cell_spec, cone_locations, varargin)
%
% Function for finding/correcting cone locations manually.  The function
% takes a datarun structure, cell specification and cone locations.  In a
% loop it plots the RFs of the rgcs given in cell_spec.  Right clicking on
% locations adds cones, left clicking removes cones, space bar (or another
% key other than "return") moves to the next rgc.  Default is for the function
% to optimize the cone location by calculating the center of mass of
% stixels around the click.  
%
% usage: new_locations = manual_cone_finding(datarun, cell_spec, cone_locations, varargin)
%
% inputs
%   datarun             standard datarun structure
%   cell_spec           specify the rgcs or types to use
%   cone_locations      Nx2 list of N cone locations
%
% optional inputs
%   gun         []            will only plot and calculate the COM on the
%                               selectec display primary
%   optimize    true    optimize the cone location that is selected by
%                       local COM
%   rad         1       radius around the clicked stixel to compute COM
%
% outputs
%   new_locations       Mx2 list of M cone locations
%   clicked_stixels     returns the stixels that were clicked
%
% GDF: 2013-02-04


% parse inputs
p = inputParser;

p.addParamValue('gun', [], @isnumeric);
p.addParamValue('optimize', true, @islogical)
p.addParamValue('rad', 1, @isnumeric);
p.addParamValue('cone_mark_color', 'c');

p.parse(varargin{:});
gun = p.Results.gun;
optimize = p.Results.optimize;
rad = p.Results.rad;
cone_mark_color = p.Results.cone_mark_color;

% FUNCTION BEGINS HERE

% get cell indices and number of cells
rgc_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(rgc_indices);

% if fits are not already present, get them -- used to zoom into RFs
if ~isfield(datarun.stas, 'fits')
    datarun = get_sta_fits_from_vision(datarun);
end

% initialize the output
new_locations = cone_locations;

% print instructions for user to display
fprintf(' left click on a stixel to add a cone \n right click on a stixel to remove a cone \n "spacebar" to move to next cell \n');

clicked_stixels = [];

clicks = 0;
% loop through RGCs
for cc = 1:num_rgcs
    
    % initialize button value
    but = 1;
    % initialize figure
    figure(1); clf;
    
    while but == 1 || but == 3

        % plot RF and cones then zoom into the receptive field
        plot_rf(datarun, datarun.cell_ids(rgc_indices(cc)), 'gun', gun)
        plot(new_locations(:,2), new_locations(:,1), 'o', 'Color', cone_mark_color);
        autozoom_to_fit(datarun, datarun.cell_ids(rgc_indices(cc)), 6, 1, 1);    

        % take input from user
        [click_xy(1), click_xy(2), but] = ginput(1);
        % round click location to nearest stixel
        click_xy = round(click_xy);

        % deal with user input: if left click, add cone,... 
        % if right click remove cone
        if but == 1
            clicks = clicks +1;
            clicked_stixels(clicks,:) = [click_xy(2), click_xy(1)];
            if optimize % optimize cone location

                % get range of values (set by rad) around the clicked point
                xrng = max(round(click_xy(2)-rad),1):min(round(click_xy(2)+rad),datarun.stimulus.field_height);
                yrng = max(round(click_xy(1)-rad),1):min(round(click_xy(1)+rad),datarun.stimulus.field_width);    

                % grab rf for cell that was clicked on
                rf = get_rf(datarun, datarun.cell_ids(rgc_indices(cc)));
                
                % only used a specified gun if called for
                if ~isempty(gun)
                    rf = rf(:,:,gun);
                else
                    rf = sum(rf,3);
                end

                % cut out an image of the cone from this cell
                cone_image = rf(xrng,yrng);

                % calculate the COM
                [new_x,new_y] = ait_centroid(cone_image);
                
                % plot some shit to see if its working
%                 if 1
%                     figure(132); clf;
%                     imagesc(cone_image); colormap gray; hold on
%                     plot(new_x, new_y, 'go')
%                     hold off
%                     pause
%                     close(gcf)
%                 end
                  
                % store the optimized COM in right new_locations
            new_locations(end+1,:) = [new_y+min(xrng)-1, new_x+min(yrng)-1];
            else
                % if not optimizing, then set cone location to that
                % clicked                
                new_locations(end+1,:) = [click_xy(2), click_xy(1)];
            end
        elseif but == 3 % handle a right click (cone removal)
            % find index of cone to be deleted
            delete_index = find(ipdm([click_xy(2), click_xy(1)],new_locations) == 0);
            if ~isempty(delete_index) 
                new_locations(delete_index,:) = [];
            else % if no cone is found report so to user and move on...
                fprintf('CLICK AGAIN: no cone found at selected point \n')
            end
        end
    end

end
