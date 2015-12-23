function the_frames = parse_frame_spec(frame_spec, num_frames)
% parse_frame_spec     Get the relevant STA frames based on the frame specification 
%
% usage:  the_frames = parse_frame_spec(frame_spec, num_frames)
%
% arguments:     frame_spec - user input
%                               ':' - return all frames
%                   positive vector - those frames (first frame is 1)
%               non-positive vector - those frames (last frame is 0)
%
%                num_frames - how many frames are in the STA
%
% outputs:     the_frames - which frames were specified (all positive numbers, first frame is 1)
%
%
%
% 2008-10 gauthier
%


if isempty(frame_spec)
    warning('empty frame specification.')
    the_frames = [];
end


switch class(frame_spec)
    case 'char'
        if strcmp(frame_spec,':')  % if :, use all frames
            the_frames = 1:num_frames;
            return
        else
            error('Frame specification ''%s'' was not recognized.',frame_spec)
        end
        
    case 'double'
        
        % convert non-positive vector to positive
        if any(frame_spec <= 0)
            frame_spec = num_frames + frame_spec;
        end

        
        % check for frames out of range

        % get list of frames which were desired but are not available
        bad_frames = setdiff(frame_spec,1:num_frames);

        % if the list is empty, no frames were bad
        if isempty(bad_frames)
            the_frames = frame_spec;
            return
            
        else % otherwise, throw an error
            error('Frame %d is not valid (must be 0 through %d).',bad_frames(1),num_frames)
        end
        
    otherwise
        fprintf('\n\n\n')
        disp(frame_spec)
        error('Above frame specification not recognized.')

end
