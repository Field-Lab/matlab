function [fittedGLM, center, STA] = glm_fit_from_block(prepped_data, grey_buffer, center, STA)

% Get the right sizes for everything
n_blocks = size(prepped_data.fitmovie,1);
block_size = size(prepped_data.fitmovie{1});
block_length = block_size(3);
grey_frames = 60;
monitor_refresh = 120;%119.5;
if grey_buffer 
    block_length = block_length+grey_frames;
    movie_mean = mean(prepped_data.fitmovie{1}(:));
    grey_movie = movie_mean*ones(block_size(1), block_size(2), grey_frames);
end % one STA length
block_size(3) = block_length*n_blocks; %movie size
fitmovie = zeros(block_size,'uint8');
fitspikes = [];

% concatenate the movie
idx = 1:block_length;
for i_block = 1:n_blocks
    if ~grey_buffer
        fitmovie(:,:,idx) = prepped_data.fitmovie{i_block};
    else
        fitmovie(:,:,idx) = cat(3, prepped_data.fitmovie{i_block}, grey_movie);
    end
    t_block_start = block_length*(i_block - 1)/monitor_refresh;
    fitspikes = [fitspikes; prepped_data.fitspikes{i_block}+t_block_start];
    idx = idx+block_length;
end

if center == 1
   [STA, center] = STA_Test(fitspikes, fitmovie, 0, 1/monitor_refresh);
end

fittedGLM = glm_fit(fitspikes, fitmovie, center, 'WN_STA', STA, 'monitor_refresh', monitor_refresh);
if isfield(prepped_data, 'testspikes')
    fittedGLM.xvalperformance = glm_predict(fittedGLM, prepped_data.testmovie, 'testspikes', prepped_data.testspikes);
end

end