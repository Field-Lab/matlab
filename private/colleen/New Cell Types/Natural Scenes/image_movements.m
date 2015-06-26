% image movements

% divide by directions in the movie
movie = 'NSEM_eye-120-3_0-3600';

cell_ids =  [337 426 486 726 736 816 919 948 950 1041 1111 1114 1115 1117 1381 2527 2617 2672 2825 2930 3081 3124 3323 3466 3634 3961 3979 4010 4027 4160 4236 4313 4564 4733 4850 4881 4988 5078 5446 5583 6218 6323 6441 6728 6756 6948 7325 7491 7564 7673 7761];
date = '2012-08-09-3';
type = 'Unknownclos';



load (['/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/', movie, '/FEM_hackreconstruct/Full_FEM.mat']);
%2012-08-09-3
% NSEM_eye-120-3_0-3600 schemeA
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';

%2012-09-27-3
% NSEM_eye-120-3_0-3600 schemeA
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';

% 2013-08-19-6
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-long-v2/fitmovie_schemeA_8pix_Identity_8pix.mat';

%2013-10-10-0
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_FEM900FF_longrast/fitmovie_schemeB_8pix_Identity_8pix.mat';

fit = fixational_EM.fitmovie;
diff_row = cell(length(fit),1);
diff_col = cell(length(fit),1);
diff = cell(length(fit),1);
diff_justangle = [];
for i = 1:length(fit)
    row = fit{i}.row_center;
    col = fit{i}.col_center;
    diff_row{i} = row - [0,row(1:end-1)];
    diff_col{i} = col - [0,col(1:end-1)];
    diff{i} = atan2(diff_col{i}, diff_row{i});
    diff{i}(diff_col{i} == 0 & diff_row{i} == 0) = nan;
    diff_justangle = [diff_justangle ; diff{i}(~isnan(diff{i}))'];
end

[n_ref, bin_ref] = histc(diff_justangle, [-pi:pi/4:3*pi/4]);



for id = 1:length(cell_ids)
    
    cell_id = cell_ids(id);
    load(['/Users/colleen/Desktop/NSEM_blocked_spikes/', date, '/BlockSpikesandRasters/organizedspikes_', type, '_', num2str(cell_id),'.mat'])
    fprintf('Computing DS STA for Cell %s... \n', num2str(cell_id))
    if id == 1
        frame_times_all_blocks = organizedspikes.block.t_frame;
        frame_times_even_blocks = frame_times_all_blocks(2:2:end);
        % frame_times = (cell2mat(frame_times_even_blocks))';
        % frame_times = frame_times(:);
        % triggers = frame_times;
        % triggers = triggers(1:100:200000);
        
        % start frame times at 0 at the beginning of the block
        %% load movie
        num_frames = 30;
        %         movie_frames = NSEMmovie.fitmovie.movie_byblock;
        
        
        % average movie
        %         movie = zeros(80,40,size(movie_frames{1}.matrix,3));
        %         for block = 1:size(movie_frames,2)
        %             movie = movie + double(movie_frames{block}.matrix)/256; % values between 0 and 1
        %         end
        %         movie_avg = mean(movie,3)'/size(movie_frames{1}.matrix,3)*num_frames;
        %
    end
    icnt = 0;
    
    spikes_by_block = organizedspikes.block.t_sp(2:2:end);
    sta=zeros(1,1,num_frames); %height, width, frames back
    sta_store = zeros(200000, 1, 30);
    for block = 1:size(spikes_by_block,1)
        frame_times = frame_times_even_blocks{block};
        
        for i=1:size(spikes_by_block{block},1)
            spike_time = spikes_by_block{block}(i);
            start=find(frame_times>spike_time,1)-num_frames;
            if(start>000)
                icnt=icnt+1;
                %             if mod(icnt, 1000) == 0
                %                 disp(icnt)
                %             end
                
                for j=1:num_frames
                    frame = double(diff{block}(start+j-1));
                    
                    %             sta(:,:,j) = sta(:,:,j) + round(reshape(F(1:3:end),width,height)'-0.5)+round(reshape(F(2:3:end),width,height)'-0.5)+round(reshape(F(3:3:end),width,height)'-0.5);
                    sta(:,:,j) = sta(:,:,j) + frame;
                    sta_store(icnt,1,j) = frame;
                    
                    
                    
                    %
                    %                 sta_store(:,:,j, icnt,1) = double(round(reshape(F(1:3:end),width,height)'-0.5));
                    %                 sta_store(:,:,j, icnt,2)=double(round(reshape(F(2:3:end),width,height)'-0.5));
                    %                 sta_store(:,:,j, icnt,3) =double(round(reshape(F(3:3:end),width,height)'-0.5));
                end
            end
        end
    end
    
    % normalize STA color
    sta = permute(sta, [2,1,3]);
    sta=sta/icnt;
    sta_store2 = sta_store(1:icnt,:,:);
    
    %     sta = sta./repmat(movie_avg,1,1,30);
    % sta = sta-repmat(movie_avg,1,1,30);
    
    sta = norm_image(sta);
    
    %     figure;
    %     for i = 1:num_frames
    %         imagesc(sta(:,:,i))
    %         colormap gray
    %         title(i)
    %         drawnow
    %         pause(0.25)
    %     end
    %
    
    %     title({date, ['Cell ', num2str(1117)], ['Frame ', num2str(i)]})
    %
    %     axis equal
    %     axis tight
    %     axis off
    temp = sta_store2(:,1,27);
    %      sta_store2(isnan(sta_store2))=[];
    DS_data = sta_store2;
    elimin = temp(isnan(temp));
    
    [n, bin] = histc(squeeze(temp(~isnan(temp))), [-pi:pi/4:3*pi/4]);
    figure;
    try
        fig = bar([-pi:pi/4:3*pi/4], n./n_ref, 'histc');
        set(gca, 'xtick', [-pi:pi/4:3*pi/4]+pi/8)
        set(gca, 'xticklabel', {'-pi','-3pi/4', '-pi/2', '-pi/4','0', 'pi/4', 'pi/2', '3pi/4'})
        title({['Cell ID: ', num2str(cell_id)],['Total Remaining Spikes: ',num2str(length(squeeze(temp(~isnan(temp)))))], [' Percent of Spikes eliminated ', num2str(length(elimin)/length(temp))]})
        dir = ['/Users/colleen/Desktop/NSEM_blocked_spikes/', date, '/NS DS/'];
        if ~exist(dir, 'dir')
            mkdir(dir);
        end
        
        save([dir, num2str(cell_id), '_NS_DS'], 'DS_data', 'n_ref', 'n')
        print(gcf,'-dpdf',sprintf('%s%s.pdf',dir,['Cell_',num2str(cell_id)]));
        
        clearvars sta_store
        clearvars sta_store2
        
    catch
        
        disp(['Probably no spikes for ' num2str(cell_id)])
    end
    
    
end



