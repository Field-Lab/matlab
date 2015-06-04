% NB 2015-05-06 
% Code for testing the input into glm_fit
% Using the same fitspikes, fitmovie, and center_coord you enter into
% glm_fit, you should be able to get an STA, with the red dot on the center
% of the cell's RF

function STA = STA_Test(fitspikes, fitmovie, center_coord)

movie_size = size(fitmovie);
STA = zeros(movie_size(1),movie_size(2),30);
fitframes = movie_size(3);

for i = 1:length(fitspikes)
    sp_frame = floor(fitspikes(i) * 120);
    if sp_frame > 29 && sp_frame<fitframes
        STA = STA+fitmovie(:,:,(sp_frame-29):sp_frame);
    end
end


for i = 1:30
   imagesc(STA(:,:,i)')
   colormap gray
   axis image
   hold on
   plot(center_coord(1), center_coord(2), 'r*')
   hold off
   pause(0.1)
end

end