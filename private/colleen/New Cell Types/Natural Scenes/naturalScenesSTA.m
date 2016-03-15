cell_ids =  [3110 6286];
date = '2013-08-19-6';
type = 'nc6';
% nc4 [322 650979 1336 1622 1641 1923 1924 1956 2182 2526 2931 7431]
% nc6 [3110 6286]
% Unknown [109 182 186 261 286 661 691 694 695 706 721 741 742 1072 1112 1156 1292 1445 1516 1668 1733 1999 2018 2105 2149 2272 2343 2372 2496 2528 2551 2751 2762 2944 3050 3196 3532 3533 3588 3635 3664 4175 4296 4308 4431 4793 4837 4925 4955 5076 5119 5135 5343 5363 5406 5446 5511 5731 6063 6256 6317 6365 6562 6591 6668 6706 6758 6781 6962 6966 6983 7007 7037 7084 7097 7309 7398 7447 7521 7563 7580 7613 7638] ]

% load(['/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-09-27-3/NSEM_mapPRJ/organizedspikes_ONPar_', num2str(cell_id),'.mat'])
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';

%2012-08-09-3
% NSEM_eye-120-3_0-3600 schemeA
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';

%2012-09-27-3
% NSEM_eye-120-3_0-3600 schemeA
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';

% 2013-08-19-6 
load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-long-v2/fitmovie_schemeA_8pix_Identity_8pix.mat';

%2013-10-10-0
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_FEM900FF_longrast/fitmovie_schemeB_8pix_Identity_8pix.mat';

for id = 1:length(cell_ids)
    
        
    cell_id = cell_ids(id);
    load(['/Users/colleen/Desktop/NSEM_blocked_spikes/', date, '/BlockSpikesandRasters/organizedspikes_', type, '_', num2str(cell_id),'.mat'])
    fprintf('Computing STA for Cell %s... \n', num2str(cell_id))
    if id == 1
        frame_times_all_blocks = organizedspikes.block.t_frame;
        frame_times_even_blocks = frame_times_all_blocks(2:2:end);
        % frame_times = (cell2mat(frame_times_even_blocks))';
        % frame_times = frame_times(:);
        % triggers = frame_times;
        % triggers = triggers(1:100:200000);
        
        % start frame times at 0 at the beginning of the block
        close 
        %% load movie
        num_frames = 30;
        movie_frames = NSEMmovie.fitmovie.movie_byblock;
        icnt = 0;
        
        % movie frames are downsampled!! FIX THIS
        % average movie
        movie = zeros(80,40,size(movie_frames{1}.matrix,3));
        for block = 1:size(movie_frames,2)
            movie = movie + double(movie_frames{block}.matrix)/256; % values between 0 and 1
        end
        movie_avg = mean(movie,3)'/size(movie_frames{1}.matrix,3)*num_frames;
        
    end

    spikes_by_block = organizedspikes.block.t_sp(2:2:end);
    sta=zeros(80,40,num_frames); %height, width, frames back
    
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
                    frame = double(movie_frames{block}.matrix(:,:,start+j))/256;
                    %             sta(:,:,j) = sta(:,:,j) + round(reshape(F(1:3:end),width,height)'-0.5)+round(reshape(F(2:3:end),width,height)'-0.5)+round(reshape(F(3:3:end),width,height)'-0.5);
                    sta(:,:,j) = sta(:,:,j) + frame;
                    
                    
                    
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
    sta = sta./repmat(movie_avg,1,1,30);
    % sta = sta-repmat(movie_avg,1,1,30);
    
    sta = norm_image(sta);
    
    figure;
    for i = 1:num_frames
        imagesc(sta(:,:,i))
        colormap gray
        title(i)
        drawnow
        pause(0.25)
    end
    
    
    save(['/Users/colleen/Desktop/NSEM_blocked_spikes/', date, '/NS STAs/', num2str(cell_id), '_NS_STA'], 'sta')
    title({date, ['Cell ', num2str(1117)], ['Frame ', num2str(i)]})
    
    axis equal
    axis tight
    axis off
end