function [inputs, refresh, duration] = get_wn_movie_ath(datarun, mdf_file, bw)

% triggers = [datarun.triggers; [datarun.triggers(end) + mean(diff(datarun.triggers)):mean(diff(datarun.triggers)):datarun.triggers(end) + 600*mean(diff(datarun.triggers))]'];
[~,height,width,duration,refresh] = get_movie_ath(mdf_file,datarun.triggers, 1,2);
mvi=load_movie(mdf_file, datarun.triggers);

if bw== 1
    inputs=zeros(height*width,duration-1);
    cnt=1;
    for j=0:duration-1
        F = round(mvi.getFrame(j).getBuffer);
        myFrames = reshape(F(1:3:end),width,height); %take out transpose
        
        to_add(:,1) = myFrames(:);
        to_add_reverse = to_add';
        inputs(:,cnt)=to_add_reverse(:);
        cnt=cnt+1;
    end
else
    
    inputs=zeros(height*width*3,duration-1);
    cnt=1;
    for j=0:duration-1
        F = mvi.getFrame(j).getBuffer;
        myFrames = reshape(F(1:3:end),width,height); %return movie with opposite orientation 
        to_add(:,1) = myFrames(:);
        
        myFrames = reshape(F(2:3:end),width,height);
        to_add(:,2) = myFrames(:);
        
        
        myFrames = reshape(F(3:3:end),width,height);
        to_add(:,3) = myFrames(:);
        
        to_add_reverse = to_add';
        inputs(:,cnt)=to_add_reverse(:);
        cnt=cnt+1;
    end
end

% inputs=(inputs*mvi.contrastValue*2)-mvi.contrastValue;
