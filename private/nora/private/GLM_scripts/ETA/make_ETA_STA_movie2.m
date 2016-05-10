threshold = 0.5;
MOVIE = 0;
STA = 1;

for i_cell=1
        eta = eta_nocone;
        sta = eta_cone;
%         minp1 = min(sta{i_cell}(:));
%         maxp1 = max(sta{i_cell}(:));
        minp = min([eta(:); eta_cone(:)]);
        maxp = max([eta(:); eta_cone(:)]);
    
    if MOVIE
        video_name = ['/Users/Nora/Desktop/RSTA_' num2str(res.cells(i_cell)) '.avi']
        writerObj = VideoWriter(video_name);
        writerObj.FrameRate = 20;
        open(writerObj);
        h = figure(1);
    end
    
    for i = 1:90
        if STA; subplot(2,1,1); end
        imagesc(eta(:,:,i)');
        axis image;
        colormap gray;
        caxis([minp maxp])
        title('RSTA')
        colorbar
        if STA
            subplot(2,1,2)
            imagesc(sta(:,:,i)');
            axis image;
            colormap gray;
            caxis([minp maxp])
            title('RSTA using cone model')
            colorbar
        end
        if MOVIE
            F = getframe(h);
            writeVideo(writerObj, F);
            clf
        else
            pause(0.1); 
        end
    end
    if MOVIE; close(writerObj); end
    
end

% save('ETAs_and_STAs_WN_rk2CP_highthres.mat', 'eta','sta')