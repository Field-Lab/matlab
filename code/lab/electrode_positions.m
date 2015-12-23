function result = electrode_positions(varargin)
% electrode_positions     return the locations of electrodes
%
% usage:
%
%   MODE 1: datarun = electrode_positions(datarun, [num_electrodes])
%
%   MODE 2: positions = electrode_positions(num_electrodes)
%
%
% arguments:     datarun - datarun struct
%         num_electrodes - how many electrodes in the array
%                           defaults to 512 if not specified in mode 1
%
% outputs:     datarun - datarun struct with positions in datarun.ei.positions
%            positions - Nx2 matrix, row E is electrode E, first column is x coordinate, second y
%
%
%  examples:
%
%  datarun = electrode_positions(datarun);
%  datarun = electrode_positions(datarun,519);
%  positions = electrode_positions(512);
%
%
%
% 2009-05  gauthier
%


% fill in default arguments
switch class(varargin{1})

    case 'struct' % mode 1
        which_mode = 1;
        datarun = varargin{1};
        switch nargin
            case 1
                num_electrodes = 512;
            case 2
                num_electrodes = varargin{2};
            otherwise
                error('too many arguments')
        end

    case 'double'  % mode 2
        which_mode = 2;
        switch nargin
            case 1
                num_electrodes = varargin{1};
            otherwise
                error('too many arguments')
        end
end


% set array_id

% if(arrayID < 500), give 64 map
% if(array >= 500 && arrayID < 1500) give 512 map
% if(array >= 1500 && < 2500) give 519 map

switch num_electrodes
    case 61
        array_id=1;
        num_positions = 64;
    case 512
        array_id=504;
        num_positions = 512;
    case 519
        array_id=1501;
        num_positions = 519;
    otherwise
        error('invalid number of electrodes (%d)',num_electrodes)
end



% load array info
array_info = load_array_info(array_id);

% get positions
positions = array_info.positions;



% return appropriately
switch which_mode
    case 1
        datarun.ei.position = positions;
        result = datarun;
    case 2
        result = positions;
end



% plot for fun
if 0
    figure
    plot(positions(:,1),positions(:,2),'.')
    axis equal
end

