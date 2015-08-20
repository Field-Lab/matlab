clear; clc
load /Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat
load /Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat
d_save =  '/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/FEM_hackreconstrunct'

%load /Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat
%d_save = '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/FEM_hackreconstruct'


if ~exist(d_save), mkdir(d_save); end

%%
% Loading stuff up

for i_blk = 1:length(NSEMmovie.fitmovie.movie_byblock
    


    block_mat = testmovie.matrix;

    testcenter_cols = [16:64];
    testcenter_rows = [16:24];
    pics = ( size(block_mat,3) ) / 120;
    i_pic = 23; i_ind = 1

    search_cols0 = [1:(80-length(testcenter_cols)+1)];
    search_rows0 = [1:(40-length(testcenter_rows)+1)];

    shifts = cell(1,pics)
    for i_pic = 1:pics
        shifted_l2.rows0  = NaN(1,120);
        shifted_l2.cols0  = NaN(1,120);
        shifted_dotsim.rows0 = NaN(1,120);
        shifted_dotsim.cols0 = NaN(1,120);


        indices = 120*(i_pic-1) + [1:120];
        reference_frame = squeeze( double(block_mat(:,:, indices(1)) )) / 255;
        imagesc(reference_frame',[0 1]); colorbar; colormap gray; %clim([0 1])

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
            display(sprintf('Pic %d, Frame %d, xval: %d, yval: %d',...
                i_pic,i_ind,  search_cols0(icol),  search_rows0(irow) ))
            shifted_l2.rows0(i_ind) = search_cols0(icol);
            shifted_l2.cols0(i_ind) = search_rows0(irow);

            max_ind = find(max_innerprod == max(max_innerprod(:)));
            [icol,irow] = ind2sub([length(search_cols0),length(search_rows0)] , max_ind);
            %display(sprintf('Pic %d, Frame %d, xval: %d, yval: %d',...
            %    i_pic,i_ind,  search_cols0(icol),  search_rows0(irow) ))
            shifted_dotsim.rows0(i_ind) = search_cols0(icol);
            shifted_dotsim.cols0(i_ind) = search_rows0(irow);


        end
        shifts{i_pic}.l2_rows = shifted_l2.rows0 -min(testcenter_rows) ;
        shifts{i_pic}.l2_cols = shifted_l2.cols0 - min(testcenter_cols);
        shifts{i_pic}.dotsim_rows = shifted_dotsim.rows0 - min(testcenter_rows);
        shifts{i_pic}.dotsim_cols = shifted_dotsim.cols0 - min(testcenter_cols);
        
        clf
        LW = 2;
        subplot(4,1,1);  hold on
        plot([1:120], shifts{i_pic}.l2_rows-16, 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('row shifts: minimize L2 error');
        ylabel('row shift from center')
        subplot(4,1,2);  hold on
        plot([1:120], shifts{i_pic}.dotsim_rows-16,'r', 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('row shifts: maximize dot similarity');
        ylabel('row shift from center')
        subplot(4,1,3);  hold on
        plot([1:120], shifts{i_pic}.l2_cols-16, 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('col shifts: minimize L2 error');
        ylabel('row shift from center')
        subplot(4,1,4);  hold on
        plot([1:120], shifts{i_pic}.dotsim_cols-16,'r', 'linewidth', LW)
        set(gca,'fontsize', 10)
        title('col shifts: maximize dot similarity');
        ylabel('rol shift from center');
        xlabel('frames (1/120) per second'); orient tall
        eval(sprintf('print %s/Block_%d_Pic_%d.pdf', d_save, i_blk,i_pic))
        
    end
end 

%{
for i_pic = 1:pics
    clf
    LW = 2;
    subplot(4,1,1);  hold on
    plot([1:120], shifts{i_pic}.l2_rows-16, 'linewidth', LW)
    set(gca,'fontsize', 10)
    title('row shifts: minimize L2 error');
    ylabel('row shift from center')
    
    subplot(4,1,2);  hold on
    plot([1:120], shifts{i_pic}.dotsim_rows-16,'r', 'linewidth', LW)
    set(gca,'fontsize', 10)
    title('row shifts: maximize dot similarity');
    ylabel('row shift from center')
    

    subplot(4,1,3);  hold on
    plot([1:120], shifts{i_pic}.l2_cols-16, 'linewidth', LW)
    set(gca,'fontsize', 10)
    title('col shifts: minimize L2 error');
    ylabel('row shift from center')
    
    subplot(4,1,4);  hold on
    plot([1:120], shifts{i_pic}.dotsim_cols-16,'r', 'linewidth', LW)
    set(gca,'fontsize', 10)
    title('col shifts: maximize dot similarity');
    ylabel('rol shift from center');
    xlabel('frames (1/120) per second'); 
    
    orient tall
    eval(sprintf('print -dpdf hello.pdf'))
end    

figure
for i_frame = 1:120
    frame = squeeze( double(block_mat(:,:,i_frame + 120*(i_pic-1)))) / 255;
    imagesc(frame',[0 1]); colormap gray; colorbar;
    pause
end
%}