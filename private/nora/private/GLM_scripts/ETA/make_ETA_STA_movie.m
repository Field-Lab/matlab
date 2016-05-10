threshold = 0.5;
STA = 0;
MOVIE = 0;

for i_cell=1:size(xval,1)
    
    disp('Making ETA')
    eta{i_cell} = RSTA(xval(i_cell,:), movie, threshold);
    
    if STA
        disp('Making STA')
        sta{i_cell} = RSTA(xval(i_cell,:), movie, 1000);
        minp = min([eta{i_cell}(:); sta{i_cell}(:)]);
        maxp = max([eta{i_cell}(:); sta{i_cell}(:)]);
    else
        minp = min(eta{i_cell}(:));
        maxp = max(eta{i_cell}(:));
    end
    
    if MOVIE
        video_name = ['/Users/Nora/Desktop/RSTA_' num2str(res.cells(i_cell)) '.avi']
        writerObj = VideoWriter(video_name);
        writerObj.FrameRate = 20;
        open(writerObj);
        h = figure(1);
    end
    
    for i = 1:90
        if STA; subplot(2,1,1); end
        imagesc(eta{i_cell}(:,:,i)');
        axis image;
        colormap gray;
        caxis([minp maxp])
        title('Residual STA')
        if STA
            subplot(2,1,2)
            imagesc(sta{i_cell}(:,:,i)');
            axis image;
            colormap gray;
            caxis([minp maxp])
            title('STA')
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