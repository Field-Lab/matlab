function center = rf_center(datarun,cell_id,center_type,coordinates)
% rf_center     retrive the center point of a RF
%
% usage:  center = rf_center(datarun,cell_id,center_type,coordinates)
%
%    last two arguments are optional
%
%
% arguments:  datarun - datarun struct
%             cell_id - cell id
%         center_type - what kind of center point to return, see options below
%         coordinates - what coordinate system to use, see options below
%
%
% outputs:     center - x,y location
%
%
% options for center_type
%
%   'sta fit'       mean of elliptical gausian fit, stored in datarun.stas.fits{cell_index}.mean
%   'com'           center of mass, stored in datarun.stas.rf_coms{cell_index}
%   'obvius'        mean of obvius 2D Gaussian fit, stored in datarun.obvius.sta_fits{cell_index}.mean
%   'vision'        mean of vision 2D Gaussian fit, stored in datarun.vision.sta_fits{cell_index}.mean
%   <function>      function applied to datarun to retreive center point
%
%       This function will always be passed two arguments: datarun, and the cell index.
%       In most cases, the function will simply specify where to retrive the center point:
%
%               @(d,cc)d.stas.special_centers{cc}
%
%       but this function allows room for unlimited complexity, including additional parameters.  for example:
%
%               @(d,cc)rf_center_complicated(d,cc,param_1,param_2)
%
%
%
%   if 'center_type' is not specified or is empty, the a default value is used
%       if 'datrun.stas.fits' exists, then 'sta_fit' is used
%       otherwise, if 'datrun.stas.rf_coms' exists, then 'com' is used
%       if neither exists, then an error is given
%
%
% options for coordinates
%
%   'sta        coordinates of the STA
%   'monitor'   coordinates of the stimulus monitor
%
%   if 'coordinates' is not specified, the coordinates of the STA is used by default
%
%
%
% gauthier  2009-04
%



% get cell index
cell_index = get_cell_indices(datarun,cell_id);


% set up defaults
if ~exist('center_type','var') || isempty(center_type)
    
    % check for datarun.stas
    if isfield(datarun,'stas')
        
        % check for datarun.stas.fits
        if isfield(datarun.stas,'fits')
            center_type = 'sta fit';
            
            % check for datarun.stas.rf_coms
        elseif isfield(datarun.stas,'rf_coms')
            center_type = 'com';
        else
            % give error
            error('datarun.stas does not have a field called ''rf_coms'' or ''fits'', so I don''t know where to find RF center points.')
        end
        
    else
        % give error
        error('datarun does not have a field ''stas'', so I don''t know where to find RF center points.')
    end
end



if ~exist('coordinates','var') || isempty(coordinates)
    coordinates = 'sta';
end


% determine how to retrieve center point
switch class(center_type)
    case 'char'
        switch center_type

            case 'com'
                % easy!
                center = datarun.stas.rf_coms{cell_index};

            case 'sta fit'
                % get sta fit mean
                center = datarun.stas.fits{cell_index}.mean;

            case 'obvius'
                % get value stored in datarun
                center = datarun.obvius.sta_fits{cell_index};

                % correct for different coordinates
                center = center + [1 1];

            case 'vision'
                % get value stored in datarun
                x = datarun.vision.sta_fits{cell_index}.mean(1);
                y = datarun.vision.sta_fits{cell_index}.mean(2);

                % correct for different coordinates
                center = [x (datarun.stimulus.field_height - y)] + 0.5;

        end

    case 'function_handle'

        % apply function
        center = center_type(datarun,cell_index);

end



% if none of the above worked out, give an error
% note: the above structure is too complicated to use an 'otherwise' case
if ~exist('center','var')
    error('center specification %s not recognized.\n',center_type)
end




% the center point is always returned in STA coordinates initially
% check whether these must be converted to monitor coordinates
switch coordinates
    case 'sta'
        % do nothing!

    case 'monitor'
        % generate transform from STA coordinates to monitor coordinates
        tform = coordinate_transform(datarun,'monitor');
        % apply it
        [xx,yy] = tformfwd(tform,center(1),center(2));
        center = [xx yy];
end


