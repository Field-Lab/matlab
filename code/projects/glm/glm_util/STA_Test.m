
% NB 2015-05-06 
function [STA, center] = STA_Test(fitspikes, fitmovie, center_verification)
%
% DESCRIPTION
% Code for testing the input into glm_fit and finding the cell's location
% 
% Using the same fitspikes and fitmovie you enter into
%   glm_fit, you should be able to get an STA using this code. 
%   The spikes need to be aligned to the triggers. 
%   See glm_fit_from_WN for an example if you don't know how
%   to do this. 
%
% INPUTS
%   fitspikes: the spike times of the neuron (in seconds)
%       These must be preperly aligned with the triggers.
%       Use STA_Test to make sure you can get an STA. Then you will know
%       that they are properly aligned.
%
%   fitmovie: the movie frame by frame. 
%       You should have a frame for every 1/120 seconds, 
%       so if the interval was two, your
%       fitmovie should have 2 of each frame
%
%   center_veritcation, logical
%       If center_verification is set the false, the code will automatically pick
%       the strongest pixel. If it is set to true, a figure will pop up asking
%       you to click on the center.
%



movie_size = size(fitmovie);
STA = zeros(movie_size(1),movie_size(2),1,30);
fitframes = movie_size(3);

for i = 1:length(fitspikes)
    sp_frame = floor(fitspikes(i) * 120);
    if sp_frame > 29 && sp_frame<fitframes
        STA = STA+reshape(double(fitmovie(:,:,(sp_frame-29):sp_frame)),[movie_size(1),movie_size(2),1,30]);
    end
end

for i = 27
   imagesc(STA(:,:,1,i)')
   colormap gray
   axis image
   title('You should see an STA here')
   pause(0.1)
end

% STA = squeeze(STA);
% STA = abs(sum(STA, 3));
% row = max(STA);
% col = max(STA');
% x = find(row == max(row));
% y = find(col == max(col));
% center = [y x];

% imagesc(STA')
% axis image
center_verification = true;
if center_verification
    title('Click on the center of the STA')
    [x, y] = ginput(1);
    center = round([x y]);
else
    hold on
    plot(center(1), center(2), 'r*')
    title('The red dot should be over the center of the STA')
    pause(0.1)
end


end
