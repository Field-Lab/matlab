% make jitter demo

stixel_size = 16;
width = 20;
true_frame = zeros(width*stixel_size, width*stixel_size, 3,10);

jitterX = load('/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Jitter/jitterX.mat');
jitterY = load('/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Jitter/jitterY.mat');


for i = 1:180*6
    x = randi([0 1],width*width,1);
    shaped_frame = reshape(x, [width,width,1]);
    
    
    sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
    
    sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
    position = [jitterX.jitterX(ceil(i/3))+1+stixel_size/2, jitterY.jitterY(ceil(i/3))+1+stixel_size/2];
    %         x and y might be reversed
    true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1), mod(i,3)+1, ceil(i/3)) = sized_frame;
    
end

for i = 1:size(true_frame,4)
    

figure; imagesc(squeeze(true_frame(:,:,:,i)))
axis equal
axis off
F(i) = getframe(gcf);
close(gcf)
end

% movie2avi(F, 'EImovie.avi')


v = VideoWriter('Jitter White Noise movie');
v.FrameRate = 15;%number_of_frames/3; 
open(v)
writeVideo(v,F);
close(v)


%% No jitter

% make jitter demo

stixel_size = 16;
width = 20;
true_frame = zeros(width*stixel_size, width*stixel_size, 3,10);

jitterX = load('/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Jitter/jitterX.mat');
jitterY = load('/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Jitter/jitterY.mat');


for i = 1:180*6
    x = randi([0 1],width*width,1);
    shaped_frame = reshape(x, [width,width,1]);
    
    
    sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
    true_frame(:,:, mod(i,3)+1, ceil(i/3)) = sized_frame;
%     sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
%     position = [1+stixel_size/2, 1+stixel_size/2];
%     %         x and y might be reversed
%     true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1), mod(i,3)+1, ceil(i/3)) = sized_frame;
%     
end

for i = 1:size(true_frame,4)
    

figure; imagesc(squeeze(true_frame(:,:,:,i)))
axis equal
axis off
F(i) = getframe(gcf);
close(gcf)
end

% movie2avi(F, 'EImovie.avi')


v = VideoWriter('No Jitter White Noise movie');
v.FrameRate = 15;%number_of_frames/3; 
open(v)
writeVideo(v,F);
close(v)