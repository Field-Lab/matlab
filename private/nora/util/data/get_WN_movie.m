function movie = get_WN_movie(mdf_file, frames)
% GET_MOVIE      Load a white noise movie

% Make fake triggers for loading
n_triggers = frames/100 + 100; % one extra for good luck
triggers = 0:(1/1.2):((1/1.2)*n_triggers);

% load movie
[mvi] = load_movie(mdf_file, triggers);

% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
% refresh = double(mvi.getRefreshTime);

if exist('frames','var')
    if frames<=duration
        tduration=frames;
    else
        error('movie too short');
    end
end
  
movie=zeros(height,width,tduration);
for i=1:tduration
    % grab movie frame
    F = mvi.getFrame(i-1).getBuffer;

    % reshape into proper image
    F = reshape(F,3,width,height);
    F = permute(F,[3 2 1]);
    movie(:,:,i) = F(:,:,1); % just take one channel since BW
end
end