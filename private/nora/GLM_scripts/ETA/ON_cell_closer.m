threshold = 0.5;

for i_cell=15
   
    
    disp('Making ETA')
    eta_CP{i_cell} = RSTA(xval_CP(i_cell,:), movie, threshold);
    disp('Making STA')
    eta_noCP{i_cell} = RSTA(xval_noCP(i_cell,:), movie, threshold);
    
    minp = min([eta_CP{i_cell}(:); eta_noCP{i_cell}(:)]);
    maxp = max([eta_CP{i_cell}(:); eta_noCP{i_cell}(:)]);
%     
%     video_name = ['/Users/Nora/Desktop/RSTA_' num2str(res.cells(i_cell)) '.avi']
%     writerObj = VideoWriter(video_name);
%     writerObj.FrameRate = 20;
%     open(writerObj);
%     h = figure(1);
    for i = 1:90
        subplot(2,1,1)
        imagesc(eta_CP{i_cell}(:,:,i)');
        axis image;
        colormap gray;
        caxis([minp maxp])
        title('RSTA with CP')
        subplot(2,1,2)
        imagesc(eta_noCP{i_cell}(:,:,i)');
        axis image;
        colormap gray;
        caxis([minp maxp])
        title('RSTA without CP')
%         F = getframe(h);
%         writeVideo(writerObj, F);
%         clf
        pause(0.05)
    end
%     close(writerObj)
    
end

% save('ETAs_and_STAs_WN_rk2CP_highthres.mat', 'eta','sta')