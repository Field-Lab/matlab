load /Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat
load /Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat 
block_mat = testmovie.matrix;


% Loading stuff up
pics = ( size(block_mat,3) ) / 120;
base_x = 36;
base_y = 16;
test_width = 8;
search_y_vals = [-15:15];
search_x_vals = [-35:35];

i_pic = 1;


% Roughly 2 seconds per second of movie  
for i_pic = 1:pics
    ems = NaN(2,120);
    indices = 120*(i_pic-1) + [1:120];
    central_pic = double(block_mat( base_x:(base_x+test_width), base_y:(base_y+test_width), indices(1) ));
    
    for i_ind = 1:120
        clear fullimage
        index = indices(i_ind);
        fullimage = double(block_mat(:,:, index));

        % ~.05 seconds per loop!
        error_onenorm = NaN(length(search_x_vals), length(search_y_vals));
        for i_x = 1:length(search_x_vals)
            for i_y = 1:length(search_y_vals)
                shift_x = search_x_vals(i_x);
                shift_y = search_y_vals(i_y);
                new_x = [ (base_x + shift_x) : (base_x + shift_x + test_width) ];
                new_y = [ (base_y + shift_y) : (base_y + shift_y + test_width) ];
                new_pic = fullimage(new_x,new_y);
                val = norm(new_pic -central_pic);
                error_onenorm(i_x,i_y) = val;
            end
        end
        
        zero_ind = find(error_onenorm == min(error_onenorm(:)));
        [ix,iy] = ind2sub([length(search_x_vals),length(search_y_vals)] , zero_ind);
        
        ems(1,i_ind) = search_x_vals(ix);
        ems(2,i_ind) = search_y_vals(iy);
        display(sprintf('Frame %d, xval: %d, yval: %d',...
            i_ind,  search_x_vals(ix),  search_y_vals(iy) ))
    end
end

for i_pic = 1:120
imagesc(squeeze(block_mat(:,:,i_pic + 2880)')); colormap gray; 
pause
end