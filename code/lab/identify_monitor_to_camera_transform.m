function T_monitor_to_camera = identify_monitor_to_camera_transform(datarun)
% identify_monitor_to_camera_transform     use datarun.piece.rig and
%               datarun.piece.optical_path_direction to determine orientation and
%               mirroring of monitor image in camera coordinates.
%
%  note: this transformation does NOT account for scaling, offset, or small rotations,
%       which must be measured using photographic mapping.
%
%
% usage:  T_monitor_to_camera = identify_monitor_to_camera_transform(datarun)
%
% arguments:     datarun - datarun struct
%
% outputs:     T_monitor_to_camera - matlab transform struct
%
%
% 2010-01  gauthier
%



% can be tested with this code:
% datarun = struct;datarun.piece.rig = 'B';datarun.piece.optical_path_direction = 'below';
% figure(1);clf;tt=[.66 .33;1 0];subplot(121);imagesc(tt);subplot(122);imagesc(imtransform(tt,identify_monitor_to_camera_transform(datarun)))


% ensure required fields exist and are nonempty
if ~isfield(datarun,'piece')
    error('information for this piece is not available.  try load_index.')
end
if ~isfield(datarun.piece,'rig') || isempty(datarun.piece.rig)
    error('Rig is not specified.  Try load_index, or set datarun.piece.rig to the correct value.')
end
if ~isfield(datarun.piece,'optical_path_direction') || isempty(datarun.piece.optical_path_direction)
    error('Optical path direction is not specified.  Try load_index, or set datarun.piece.optical_path_direction to the correct value (should be "above" or "below").')
end


switch [datarun.piece.rig ' ' datarun.piece.optical_path_direction]
    case 'A above'
        T_monitor_to_camera = maketform('affine',[0 1;0 0;1 0],[-1 0;0 0;0 1]);
    case 'A below'
        T_monitor_to_camera = maketform('affine',[0 1;0 0;1 0],[1 0;0 0;0 1]);
    case 'B above'
        T_monitor_to_camera = maketform('affine',[0 1;0 0;1 0],[1 0;0 0;0 -1]);
    case 'B below'
        T_monitor_to_camera = maketform('affine',[0 1;0 0;1 0],[-1 0;0 0;0 -1]);
    case 'C above'
        T_monitor_to_camera = maketform('affine',[0 1;0 0;1 0],[-1 0;0 0;0 -1]);
    otherwise
        error('rig and optical path ''%s'' not recognized',[datarun.piece.rig ' ' datarun.piece.optical_path_direction])
end

