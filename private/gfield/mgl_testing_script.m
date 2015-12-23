%% display a white noise movie
screen = 2;
display_params = mglResolution(screen);

% Initialize display and set to zero.
mglOpen(screen);

mglClearScreen([0.5 0.5 0.5])
mglScreenCoordinates
mglFlush;
mglClearScreen([0.5 0.5 0.5])
mglFlush;


% setup screen parameters
stixel_width = 3; stixel_height = 3; % number of pixels per stixel
num_width = 200; num_height = 200; % number of stixels in x and y
stim_width = stixel_width * num_width;
stim_height = stixel_height * num_height;
x_begin = display_params.screenWidth/2 - stim_width/2; % location of screen on display
y_begin = display_params.screenHeight/2 - stim_height /2;
max_run_time = 5; % duration the stimulus will run in seconds

interval = 1; % number of monitor frames per stimulus frame
RGB = true; % true = make RGB stimulus, false = make BW stimulus
preformat = true; % preformat the noise_matrix to be in uint8 rgba x imageWidth x imageHeight

% array that keeps the time it takes to run various things
profileTimes = nan(4,max_run_time*mglGetParam('frameRate'));
profileNames = {'randi','mglCreateTexture','mglBltTexture','mglFlush etc'};

% initialize large matrices
%base_matrix = zeros(num_width,num_width,3);
%noise_matrix = matrix_scaled_up(base_matrix, stixel_width);

startTime = mglGetSecs;
frame_counter = 1;
while mglGetSecs(startTime) < max_run_time;
    if frame_counter * interval/display_params.frameRate <= mglGetSecs(startTime)
        
        profileStartTime = mglGetSecs;
        if RGB % rgb white noise
            base_matrix = randi([0,1], num_height, num_width, 3)*255;
            noise_matrix = base_matrix;
            noise_matrix(:,:,4) = 255;
        end
        
        profileTimes(1,frame_counter) = mglGetSecs(profileStartTime);
        profileStartTime = mglGetSecs;
        
	% mglCreateTexture
        texture = mglCreateTexture(noise_matrix,[],0,{'GL_TEXTURE_MAG_FILTER','GL_NEAREST'});
    	profileTimes(2,frame_counter) = mglGetSecs(profileStartTime);
    	profileStartTime = mglGetSecs;
	
	% mglBltTexture
        mglBltTexture(texture, [x_begin, y_begin, stim_width, stim_height], -1, -1, 0);
	profileTimes(3,frame_counter) = mglGetSecs(profileStartTime);
	profileStartTime = mglGetSecs;
	
	% mglFlush
        mglFlush;
        mglDeleteTexture(texture);

        % to abort stimulus presentation
        if mglGetKeys([66,65]) % keys "a" and "b"
            mglClearScreen(0.5)
            mglFlush;
            break
        end
	profileTimes(4,frame_counter) = mglGetSecs(profileStartTime);
        frame_counter = frame_counter + 1;
	
    end
end

if preformat
  disp(sprintf('%ix%i base %ix%i stixel (uint8 precomputed matrix):',num_width,num_height,stixel_width,stixel_height));
else
  disp(sprintf('%ix%i base %ix%i stixel:',num_width,num_height,stixel_width,stixel_height));
end
% display the median time 
for i = 1:size(profileTimes,1)
  disp(sprintf('%s: %0.2fms',profileNames{i},1000*median(profileTimes(i,1:frame_counter-1))));
end

% clear screen and print to command window the number of frames displayed
mglClearScreen(0.5)
mglFlush;
disp(sprintf('Number of frames: %i', frame_counter));