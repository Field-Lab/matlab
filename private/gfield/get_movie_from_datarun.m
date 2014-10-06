function [mov,height,width,duration,refresh] = get_movie_from_datarun(datarun, frames)
% GET_MOVIE      Load a white noise movie
%
%   length(triggers) must be >= 10*stimulus_interval(1,2...)
%
%greschner
%
% modified by gdf to first check whether datarun.stimulus.java_movie is
% filled. If it is, then the movie object is taken from there.
% gdf 2011-12-02

mvi = datarun.stimulus.java_movie;

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
  
mov=zeros(height,width,3,tduration);
for i=1:tduration
    % grab movie frame
    F = mvi.getFrame(i-1).getBuffer;

    % reshape into proper image
    F = reshape(F,3,width,height);
    mov(:,:,:,i) = permute(F,[3 2 1]);
end


%test:
%imagesc(mov(:,:,1,1)) same orientation as on stimulus monitor

