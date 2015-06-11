clear; clc

runoptions.compute_traces     = false;
runoptions.organize_singlemat = true;

%movie_name = 'NSEM_eye-120-3_0-3600';
movie_name = 'NSEM_eye-long-v2';
%movie_name = 'NSEM_FEM900FF_longrast';
stimdir    = '/Users/akheitman/NSEM_Home/Stimuli';
%stimdir    = '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli';

scheme = 'schemeA';
if strcmp(movie_name, 'NSEM_FEM900FF_longrast')
    scheme = 'schemeB';
end

eval(sprintf('load %s/%s/testmovie_%s_8pix_Identity_8pix.mat',...
    stimdir, movie_name,scheme));

if runoptions.compute_traces
eval(sprintf('load %s/%s/fitmovie_%s_8pix_Identity_8pix.mat',...
    stimdir, movie_name,scheme));
end

d_save = sprintf('%s/%s/FEM_hackreconstruct', stimdir, movie_name')
if ~exist(d_save), mkdir(d_save); end
image_dir = sprintf('%s/traces', d_save);
if ~exist(image_dir), mkdir(image_dir); end

%%
if runoptions.compute_traces
for i_blkind = 0:length(NSEMmovie.fitmovie.movie_byblock)
    
    if i_blkind == 0
         stimblock = 1;
         block_mat = testmovie.matrix;
    else
        stimblock = NSEMmovie.fitmovie.ind_to_block(i_blkind);
        block_mat = NSEMmovie.fitmovie.movie_byblock{i_blkind}.matrix;
    end
    
    testcenter_cols = [16:64];
    testcenter_rows = [16:24];
    pics = ( size(block_mat,3) ) / 120;
    search_cols0 = [1:(80-length(testcenter_cols)+1)];
    search_rows0 = [1:(40-length(testcenter_rows)+1)];
    clear fem
    if stimblock == 1
        fem.movietype = 'testmovie';
        fem.stimblock = stimblock;
    else
        fem.movietype = 'fitmovie';
        fem.stimblock = stimblock;
        fem.fitmovie_index  = i_blkind;
    end
    fem.shifts = cell(1,pics);
    for i_pic = 1:pics
        display(sprintf('Working on Block %d, Image %d',  stimblock,i_pic))
        shifted_l2.rows0  = NaN(1,120);
        shifted_l2.cols0  = NaN(1,120);
        shifted_dotsim.rows0 = NaN(1,120);
        shifted_dotsim.cols0 = NaN(1,120);
        indices = 120*(i_pic-1) + [1:120];
        reference_frame = squeeze( double(block_mat(:,:, indices(1)) )) / 255;
       % imagesc(reference_frame',[0 1]); colorbar; colormap gray; %clim([0 1])

        for i_ind = 1:120
            clear fullimage
            index = indices(i_ind);
            newframe = squeeze( double(block_mat(:,:, index) )) / 255;
            newframe_center = newframe(testcenter_cols,testcenter_rows);

            % ~.05 seconds per loop!
            error_l2norm  = NaN(length(search_cols0), length(search_rows0));
            max_innerprod = NaN(length(search_cols0), length(search_rows0));
            for i_col = 1:length(search_cols0)
                for i_row = 1:length(search_rows0)
                    ref_col = [search_cols0(i_col) :(search_cols0(i_col)+length(testcenter_cols)-1)];
                    ref_row = [search_rows0(i_row) :(search_rows0(i_row)+length(testcenter_rows)-1)];
                    ref_newcenter = reference_frame(ref_col,ref_row);
                    error_l2norm(i_col,i_row) = norm(newframe_center-ref_newcenter);
                    newdot = sum(sum(ref_newcenter .* newframe_center) );
                    max_innerprod(i_col,i_row) = newdot / ( norm(ref_newcenter) * norm(newframe_center) );
                end
            end

            zero_ind = find(error_l2norm == min(error_l2norm(:)));
            [icol,irow] = ind2sub([length(search_cols0),length(search_rows0)] , zero_ind);
            
            shifted_l2.rows0(i_ind) = search_cols0(icol);
            shifted_l2.cols0(i_ind) = search_rows0(irow);

            max_ind = find(max_innerprod == max(max_innerprod(:)));
            [icol,irow] = ind2sub([length(search_cols0),length(search_rows0)] , max_ind);
            shifted_dotsim.rows0(i_ind) = search_cols0(icol);
            shifted_dotsim.cols0(i_ind) = search_rows0(irow);

            %display(sprintf('Pic %d, Frame %d, xval: %d, yval: %d',...
             %   i_pic,i_ind,  search_cols0(icol),  search_rows0(irow) ))
        end
        fem.shifts{i_pic}.l2_rows = shifted_l2.rows0 -min(testcenter_rows) ;
        fem.shifts{i_pic}.l2_cols = shifted_l2.cols0 - min(testcenter_cols);
        fem.shifts{i_pic}.dotsim_rows = shifted_dotsim.rows0 - min(testcenter_rows);
        fem.shifts{i_pic}.dotsim_cols = shifted_dotsim.cols0 - min(testcenter_cols);
        
        clf
        subplot(3,1,1); 
        imagesc(reference_frame',[0 1]); colorbar; colormap gray; hold on; %clim([0 1])
        title(sprintf('Block: %d, Image: %d',stimblock,i_pic),'interpreter','none')
        
        LW = 2;
        subplot(6,1,3);  hold on
        plot([1:120], fem.shifts{i_pic}.l2_rows, 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('row shifts: minimize L2 error');
        ylabel('pixels')
        subplot(6,1,4);  hold on
        plot([1:120], fem.shifts{i_pic}.dotsim_rows,'r', 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('row shifts: maximize dot similarity');
        ylabel('pixels')
        subplot(6,1,5);  hold on
        plot([1:120], fem.shifts{i_pic}.l2_cols, 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('col shifts: minimize L2 error');
        ylabel('pixels');
        subplot(6,1,6);  hold on
        plot([1:120], fem.shifts{i_pic}.dotsim_cols,'r', 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('col shifts: maximize dot similarity');
        ylabel('pixels');
        xlabel('frames (1/120) per second'); 
        
        if stimblock == 1
            orient tall
            eval(sprintf('print -dpdf %s/TestMovie_Image_%d.pdf', image_dir,i_pic))
        else
            orient tall
            eval(sprintf('print -dpdf %s/FitMovie_Block_%d_Image_%d.pdf', image_dir,stimblock,i_pic))
        end
    
    end
    if stimblock == 1
        eval(sprintf('save %s/fem_testmovie.mat fem',  d_save))
    else
        eval(sprintf('save %s/fem_block_%d.mat fem',  d_save, stimblock))
    end
end 
end
%%

if runoptions.organize_singlemat
 
fitblocks  = testmovie.params.fitblocks;
testblocks = testmovie.params.rasterblocks;


fixational_EM.moviename = movie_name;
fixational_EM.note1   = 'reconstructed eye-movements from downsampled stim';
fixational_EM.note2   = 'AKHeitman -- fairly confident this is correct';
fixational_EM.computetime = datestr(clock);
fixational_EM.computemfile = mfilename;
fixational_EM.fitmovie = cell(1,length(fitblocks))
for i_blkind = 0:length(fitblocks)
    if i_blkind == 0
        stimblock = 1;
        eval(sprintf('load %s/fem_testmovie.mat fem',  d_save))
    else
        stimblock = fitblocks(i_blkind);
        eval(sprintf('load %s/fem_block_%d.mat fem',  d_save, stimblock))
    end
    
    images = length(fem.shifts);
    fpi = 120;  % 120 frames per image
    row_center   = NaN(1, images*fpi);
    col_center   = NaN(1, images*fpi);
    displacement = NaN(1,images*fpi);    
    for i_image = 1:images
        indices  = (i_image-1)*fpi + [1:fpi];
        new_rows = .5* (fem.shifts{i_image}.l2_rows + fem.shifts{i_image}.dotsim_rows);
        new_cols = .5* (fem.shifts{i_image}.l2_cols + fem.shifts{i_image}.dotsim_cols);
        
        change_rows = [0 diff(new_rows)];
        change_cols = [0 diff(new_cols)];
        change_rad  = ([change_rows.^2 + change_cols.^2]).^(.5);
        
        row_center(indices) = new_rows;
        col_center(indices) = new_cols;
        displacement(indices) = change_rad;
    end
    clear dummy
    dummy.stimblock  = stimblock;
    dummy.row_center = row_center;
    dummy.col_center = col_center;
    dummy.displacement = displacement;
    dummy.displacement_note = 'net downsamled pixel change from previous frame';
    if i_blkind == 0
        dummy.stimblock = fitblocks;
        fixational_EM.testmovie = dummy;
    else
        fixational_EM.fitmovie{i_blkind} = dummy;
    end
    display(sprintf('Done with block %d', i_blkind));
end

eval(sprintf('save %s/Full_FEM.mat fixational_EM', d_save))


end
