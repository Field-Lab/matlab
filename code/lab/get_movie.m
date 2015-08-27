function [mov,height,width,duration,refresh] = get_movie(mdf_file, triggers, frames)
% GET_MOVIE      Load a white noise movie
%
%   length(triggers) must be >= 10*stimulus_interval(1,2...)
%
%greschner

% load movie
[mvi] = load_movie(mdf_file, triggers);

% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);

if exist('frames','var')
    if frames<=duration
        tduration=frames;
    else
        error('movie too short');
    end
end
  

if length(tduration) > 1
    mov=zeros(height,width,3,length(tduration));
    cntr = 0;
    for i=tduration(1):tduration(length(tduration))
        F = mvi.getFrame(i-1).getBuffer;
        
        F = reshape(F,3,width,height);
        cntr = cntr + 1;
        mov(:,:,:,cntr) = permute(F,[3 2 1]);
    end
else
    mov=zeros(height,width,3,tduration);
    for i=1:tduration
        % grab movie frame
        F = mvi.getFrame(i-1).getBuffer;

        % reshape into proper image
        F = reshape(F,3,width,height);
        mov(:,:,:,i) = permute(F,[3 2 1]);
    end
end


%test:
%imagesc(mov(:,:,1,1)) same orientation as on stimulus monitor

