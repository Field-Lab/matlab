function [inputs, refresh, duration] = get_wn_movie_ath_rgb(datarun, movie_name)

triggers=datarun.triggers; % onsets of the stimulus presentation in s
mdf_file=['/Volumes/Analysis/stimuli/white-noise-xml/', movie_name];
[~,height,width,duration,refresh] = get_movie_ath(mdf_file,triggers, 1,2);
mvi=load_movie(mdf_file, triggers);

inputs=zeros(height*width,3,duration-1);
cnt=1;
for j=0:duration-1
    F = round(mvi.getFrame(j).getBuffer);
    myFrames = reshape(F(1:3:end),width,height)';    
    inputs(:,1,cnt)=myFrames(:);
    myFrames = reshape(F(2:3:end),width,height)';
    inputs(:,2,cnt)=myFrames(:);
    myFrames = reshape(F(3:3:end),width,height)';
    inputs(:,3,cnt)=myFrames(:);
    cnt=cnt+1;
end
inputs=(inputs*0.96)-0.48;
