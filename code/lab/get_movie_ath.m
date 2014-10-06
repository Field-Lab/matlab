function [mov,height,width,duration,refresh] = get_movie_ath(mdf_file, triggers, start_frame, end_frame)
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

if exist('end_frame','var')
    if end_frame<=duration
        tduration=end_frame;
    else
        error('movie too short');
    end
end
  
mov=zeros(height,width,3,tduration);
for i=start_frame:tduration
    % grab movie frame
    F = mvi.getFrame(i-1).getBuffer;

    % reshape into proper image
    F = reshape(F,3,width,height);
    mov(:,:,:,i-start_frame+1) = permute(F,[3 2 1]);
end

%test:
%imagesc(mov(:,:,1,1)) same orientation as on stimulus monitor

