sdfdsjhkbhjkglk% Code to make stimulus to just have white noise over the large cells and a
% gray background everywhere else
file_path = '/Users/colleen/matlab/private/colleen/New Cell Types/Stimulus Code/test.m';
screen_size_y = 320; % vertical size
screen_size_x = 320; % hortizontal size
stixel_size = 20;
large_cell_stixels = [10,5; 10,10; 2,12]; % horizontal over from top left then vertical down from top corner
if max(large_cell_stixels(:,1)*stixel_size) > screen_size_x
    disp('error: x dimension of chosen stixels doesn''t fit on screen')
    return
elseif max(large_cell_stixels(:,2)*stixel_size) > screen_size_y
    disp('error: y dimension of chosen stixels doesn''t fit on screen')
    return
end

    

myMap = zeros(screen_size_y, screen_size_x); % pixesl on the screen
scale_factor_x = screen_size_x/ stixel_size;
scale_factor_y = screen_size_y/stixel_size;

for i = 1:size(large_cell_stixels,1)
    pix_y = (stixel_size*large_cell_stixels(i, 1)-(stixel_size-1)):(stixel_size)*large_cell_stixels(i,1);
    pix_x = (stixel_size*large_cell_stixels(i,2)-(stixel_size-1)):(stixel_size)*large_cell_stixels(i,2);
    myMap(pix_x,pix_y) = i;
end



dlmwrite(file_path, myMap, 'delimiter', '\t', 'newline', 'pc');
savedMap = dlmread(file_path);
figure
imagesc(savedMap)
axis equal