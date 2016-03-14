function [inputs, refresh, duration] = get_wn_movie_ath(datarun, movie_name)

triggers=datarun.triggers; % onsets of the stimulus presentation in s
mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', movie_name];

mvi=load_movie(mdf_file, triggers);
% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);
% 
% if exist('end_frame','var')
%     if end_frame<=duration
%         tduration=end_frame;
%     else
%         error('movie too short');
%     end
% end
  

% inputs=zeros(height*width,duration-1);
% cnt=1;
% for j=0:duration-1
%     F = round(mvi.getFrame(j).getBuffer);
%     myFrames = reshape(F(1:3:end),width,height)';    
%     inputs(:,cnt)=myFrames(:);
%     cnt=cnt+1;
% end
% inputs=(inputs*0.96)-0.48;

inputs=zeros(height,width,duration);
for j=0:duration-1
    if mod(j, 1000)==0
        j
    end
    F = round(mvi.getFrame(j).getBuffer);
    myFrames = reshape(F(1:3:end),width,height);    
    myFrames(myFrames==0) = -.48;
    myFrames(myFrames==1) = .48;
    inputs(:, :, j+1)=myFrames';
end
