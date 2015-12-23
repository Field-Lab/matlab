function varargout = load_array_info(varargin)
% array_info     get information about an array
%
%
% arguments:    source - EITHER the array id OR datarun struct with field datarun.piece.array_id
%
%
%
%   MODE 1: datarun = load_array_info(datarun)
%
% outputs:    datarun - datarun struct with fields
%               datarun.ei.position - Nx2 matrix of electrode positions, each row an electrode
%             array - struct with fields; see below
%
%
%
%   MODE 2: array_info = load_array_info(< datarun | array_id >)
%
% outputs:     array - struct with fields
%                         positions    - Nx2 matrix, electrode positions
%                         disconnected - N-length binary vector, true = disconnected
%                         corners      - Mx2 matrix, positions of array "corners"  NOTE this is not the same as corner 
%                                           electrode positions for blunt corners!
%                         x_bounds     - 2-length vector, maximum electrode extent, plus a little padding to make a pretty plot
%                         y_bounds     -  "
%                         image        - matrix, image of the array (screen capture from platinization)
%                                           the connector is on the left, and the array is seen from above
%                         T_array_to_array_image -
%                                        transformation from array coordinates to array image coordinates.
%                                        used to more easily selected alignment points between array and camera image.
%                         T_array_to_vision_ei - 
%                                        transformation from array coordinates to vision EI coordinates
%                         shape        - string, shape of array ('rectangle' or 'hexagon')
%                         spacing      - scalar, electrode spacing in microns (30, 60, or 120)
%                         electrodeMap - java object of the electrode map, useful for finding neighboring electrodes
%              array - returns the same thing twice!  This is kludgy, but
%                      makes the second return always consistent. 
%
%
% NOTE: can also specify mode with a second argument, e.g.
%
%   array_info = load_array_info(datarun,2)
%
%
%
%  The output of this function can be tested with the script test_load_array_info.m
%
%
%
% 2010-01  gauthier
% 2010-05  phli - modified to pull array_id from globals file if necessary
% 2010-05  phli - modified so that will return full array info as second arg in mode 1
% 2010-09  phli - abstracted out array image loading; this is independently useful
%



% GET ARRAY ID AND IDENTIFY MODE


% if a struct was provided...
if isstruct(varargin{1})

    % check for appropriate field...
    if isfield(varargin{1},'piece') && isfield(varargin{1}.piece,'array_id') && ~isempty(varargin{1}.piece.array_id)
        array_id = varargin{1}.piece.array_id;
    else
        datarun = load_globals(varargin{1});
        if isfield(datarun, 'globals') && datarun.globals.getImageCalibrationParams.arrayID
            array_id = datarun.globals.getImageCalibrationParams.arrayID;
            fprintf('!!!!!! ARRAY ID WAS LOADED FROM GLOBALS: %d !!!!!!', array_id)
        else
            error('datarun struct does not have field ''datarun.piece.array_id''.')
        end
    end

    % default to returning datarun struct
    which_mode = 1;


elseif isnumeric(varargin{1}) % if a number was provided...

    % assume it's the array ID
    array_id = varargin{1};

    % default to returning array info struct
    which_mode = 2;


else
    % if input is neither struct nor numeric, give error
    disp(varargin{1})
    error('array specification not recognized (printed above)')
end



% if there's a second argument...
if nargin == 2  

    % it will specify the mode, overriding the default
    if varargin{2} == 1 || varargin{2} == 2
        which_mode = varargin{2};
    else
        error('mode type must be 1 or 2')
    end

elseif nargin > 2 % shame the fool that provides more than 2 arguments
    error('can not accept > 2 arguments')
end



% handle special case: disallow mode 1 if the datarun struct was not provided
if ~isstruct(varargin{1}) && which_mode == 1
    error('in mode 1, user must provide datarun struct')
end



% GET ELECTRODE POSITIONS

% get electrodeMap object
electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(array_id);


% get location of each electrode
num_positions = electrodeMap.getNumberOfElectrodes - 1; % Don't include the first trigger electrode
positions = zeros(num_positions,2);
for pp = 1:num_positions
    positions(pp,1) = electrodeMap.getXPosition(pp);
    positions(pp,2) = electrodeMap.getYPosition(pp);
end





% USE ARRAY ID TO IDENTIFY ARRAY PARAMETERS

% if(arrayID < 500), give 64 map
% if(array >= 500 && arrayID < 1500) give 512 map
% if(array >= 1500 && < 2500) give 519 map

if array_id < 500  % 64 electrode array
    corner_electrodes = [8 19 30 40 51 62];
    x_bounds = 350*[-1 1];
    y_bounds = 300*[-1 1];
    shape = 'hexagon';
    spacing = 60;
    if vision_version >= 7002000
        T_array_to_array_image = maketform('affine', [-120 -240; -240 0; -120 240], [217.5 253.5; 337.5 193.5; 457.5 253.5]);
    else
        T_array_to_array_image = maketform('affine', [120 -240; 240 0; 120 240], [217.5 253.5; 337.5 193.5; 457.5 253.5]);
    end
    T_array_to_vision_ei = maketform('affine',[0 1;0 0;1 0],[0 -1;0 0;1 0]); % vertical flip

elseif (array_id >= 500) && (array_id < 1500) % 512 array
    corners = [-945 450;945 450;945 -450;-945 -450];
    corner_electrodes = [249 392 512 129];
    x_bounds = 1100*[-1 1];
    y_bounds = 500*[-1 1];
    shape = 'rectangle';
    spacing = 60;
    %array_image = imrotate(array_image,270,'bilinear');  % 2010-08-16; replaced image so now shouldn't need this rotation % 2010-10-07 Now the array_image logic is abstracted to somewhere else, so this line is really out of date, but keep for now...
    T_array_to_array_image = maketform('affine',[-1288.75 -951.25; -1288.75 966.25 ; 1276.6552734375 966.25 ],[1 1;768 1;768 1027]);
    T_array_to_vision_ei = maketform('affine',[0 1;0 0;1 0],[0 -1;0 0;-1 0]); % 180 deg rotation

elseif array_id == 1530 || array_id == 1501 || array_id == 1504 % 519, 30 �m array

%    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
%    fprintf('!!!!!!!!   WARNING: array IDs have not been set yet for values > 1500   !!!!!!!!!!!!\n')
%    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')

    corner_electrodes = [264 195 126 4 455 386];
    x_bounds = 900*[-1 1];
    y_bounds = 800*[-1 1];
    shape = 'hexagon';
    spacing = 30;

    if vision_version >= 7002000
        T_array_to_array_image = cp2tform([-180 105; 0 195; 180 105; 180 -105; 0 -195; -180 -105], [658.5, 121; 803 409; 658.5 697; 322.5 697; 179 409; 322.5 121], 'affine');
    else
        %T_array_to_array_image = maketform('affine',[210 -360; -210 -360; 210 360],[664.5 702.5; 328.5 702.5; 664.5 126.5]);
        T_array_to_array_image = maketform('affine',[-210 360; 210 360; -390 0 ],[322.5 120.5; 658.5 120.5; 178.5 408.5]);
    end

    T_array_to_vision_ei = maketform('affine',[0 1;0 0;1 0],[0 -1;0 0;-1 0]); % 180 deg rotation

elseif array_id == 1512  % 519, 120 �m array

%    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
%    fprintf('!!!!!!!!   WARNING: array IDs have not been set yet for values > 1500   !!!!!!!!!!!!\n')
%    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')

    corner_electrodes = [264 195 125 4 455 386];
    x_bounds = 900*[-1 1];
    y_bounds = 800*[-1 1];
    shape = 'hexagon';
    spacing = 120;

    if vision_version >= 7002000
        T_array_to_array_image = cp2tform([-720 420; 0 780; 720 420; 720 -420; 0 -780; -720 -420], [658.5, 121; 803 409; 658.5 697; 322.5 697; 179 409; 322.5 121], 'affine');
    else
        T_array_to_array_image = maketform('affine',[840 -1440; -840 -1440; 840 1440],[658.5 696.5; 322.5 696.5; 658.5 120.5]);
    end
    
    T_array_to_vision_ei = maketform('affine',[0 1;0 0;1 0],[0 -1;0 0;-1 0]); % 180 deg rotation

elseif (array_id >= 1500) && (array_id < 2500) % 519 array

    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
    fprintf('!!!!!!!!   WARNING: array IDs have not been set yet for values >= 1500   !!!!!!!!!!!\n')
    fprintf('!!!!!!!!   In addition, this is an unknown > 1500 array id;              !!!!!!!!!!!\n')
    fprintf('!!!!!!!!       transforms may be broken                                  !!!!!!!!!!!\n')
    fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')

    num_positions = 519;
    corner_electrodes = [264 195 125 4 455 386];
    x_bounds = 900*[-1 1];
    y_bounds = 800*[-1 1];
    shape = 'hexagon';
    spacing = 60;
    T_array_to_array_image = maketform('affine',[0 1;0 0;1 0],[0 1;0 0;1 0]); % identity transformation -- WRONG!
    T_array_to_vision_ei = maketform('affine',[0 1;0 0;1 0],[0 -1;0 0;-1 0]); % 180 deg rotation

else
    error('array id  %d not recognized.',array_id)
end


% Get array image
array_image_file = fullfile(matlab_code_path(), get_array_image(array_id));
array_image = imread(array_image_file);


% IDENTIFY DISCONNECTED ELECTRODES

% initialize storage variable
disconnected = false(size(positions,1),1);

% for each electrode
for ee=1:size(positions,1)
    % determine whether disconnected
    if electrodeMap.isDisconnected(ee)==1
        disconnected(ee) = true;
    end
end





% SCALE IF NEEDED
if vision_version < 7002000 % Pitch has been corrected from 7.2.0 forward
    if array_id == 1530 || array_id == 1501 % 519, 30 �m array
        positions = positions / 2;  % IS THIS RIGHT????
        x_bounds = x_bounds / 2;
        y_bounds = y_bounds / 2;
        
    elseif array_id == 1512  % 519, 120 �m array
        positions = positions * 2;  % IS THIS RIGHT????
        x_bounds = x_bounds * 2;
        y_bounds = y_bounds * 2;
    end
elseif array_id == 1512 % But this is still fucked up because 1512 is not supposed to be a 120 um...
    positions = positions .* 4;
    x_bounds = x_bounds .* 4;
    y_bounds = y_bounds .* 4;
end


% GET CORNERS
if ~exist('corners','var')
    corners = positions(corner_electrodes,:);
end



% plot for fun
if 0
    % plot electrode points
    figure(3);
    %imagesc(array_image);axis image
    hold on
    plot(positions(:,1),positions(:,2),'.')
    hold on
    axis equal
    % plot electrode IDs
    for ee=1:size(positions,1)
        text(positions(ee,1),positions(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',12,...
            'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
    % plot outline
    temp = [corners;corners(1,:)];
    plot(temp(:,1),temp(:,2),'r')
end




% load stuff into struct
array_info = struct;
array_info.array_id = array_id;
array_info.positions = positions;
array_info.disconnected = disconnected;
array_info.x_bounds = x_bounds;
array_info.y_bounds = y_bounds;
array_info.corners = corners;
array_info.corner_electrodes = corner_electrodes;
array_info.image_file = array_image_file;
array_info.image = array_image;
array_info.T_array_to_array_image = T_array_to_array_image;
array_info.T_array_to_vision_ei = T_array_to_vision_ei;
array_info.shape = shape;
array_info.spacing = spacing;
array_info.electrodeMap = electrodeMap;

% return appropriately
% 
switch which_mode
    case 1
        datarun = varargin{1};
        datarun.ei.position = positions;
        datarun.ei.array_bounds_x = x_bounds;
        datarun.ei.array_bounds_y = y_bounds;
        varargout{1} = datarun;
        varargout{2} = array_info;
    case 2
        varargout{1} = array_info;
        varargout{2} = array_info;
end
