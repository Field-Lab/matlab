%% Add paths

% Change this to your github repo
repo_location = 'Documents/MATLAB'; 

% Then these will work
javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar')
addpath(genpath([repo_location '/matlab/code/lab']))
addpath(genpath([repo_location '/matlab/code/projects/glm']))

%% Fit the GLM to a white noise run
alex_test = glm_fit_from_WN(800, '2015-03-09-2/d05-27-norefit/data008-from-d05-d27/data008-from-d05-d27', 'BW-10-8-0.48-11111-32x32');

%% Look at the model output
plotfilters(alex_test)

%% Load up the Natural Scenes test movie
idx = 1:120;
for i = 1:30
   load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_' num2str(i) '.mat']); 
   testmovie(:,:,idx) = movie;
   idx = idx+120;
end
% Downsample and take the middle to match the white noise run
testmovie_cut = imresize(testmovie(81:240,:,:), 0.2); 

%% Make predictions!
xval = glm_predict(alex_test, testmovie_cut);
plotrasters(xval, alex_test)