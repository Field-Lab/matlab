number_of_frames = 10*30;
for i = 1:number_of_frames
% rand_int = round(rand(16,16,3));
rand_int = normrnd(0.5,0.16,[16,16, 3]);


figure; imagesc(rand_int)
axis equal
axis off
F(i) = getframe(gcf);
close(gcf)
end

% movie2avi(F, 'EImovie.avi')


v = VideoWriter('Gaussian White Noise movie');
v.FrameRate = 30;%number_of_frames/3; 
open(v)
writeVideo(v,F);
close(v)

