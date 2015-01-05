cell=2;

hold on
plot(datarun.vision.timecourses(cell).r,'r')
plot(datarun.vision.timecourses(cell).g,'g')
plot(datarun.vision.timecourses(cell).b,'b')

RGB_linear=[1, linsolve(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).g),...
linsolve(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).b)];

%%
hold on
scatter(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).g,'g')
scatter(datarun.vision.timecourses(cell).r,datarun.vision.timecourses(cell).b,'b')
vec=[min(datarun.vision.timecourses(cell).r) max(datarun.vision.timecourses(cell).r)];
plot(vec, RGB_linear(2)*vec,'g');
plot(vec, RGB_linear(3)*vec,'b');

%% 
RGB=[datarun.vision.timecourses(cell).r ...
    datarun.vision.timecourses(cell).g ...
    datarun.vision.timecourses(cell).b];

[U,S,V]=svd(RGB);
RGB_SVD=V(:,1);

%% Spikes times in terms of frames
t_frame = t_frame_interpAH(datarun.triggers);
spikes=datarun.spikes{2};
spikes_frame=zeros(size(spikes));
for i=1:length(spikes)
    spikes_frame(i)=find(t_frame==max(t_frame(t_frame<spikes(i))));
end

%% How I do it in the code

spikes=datarun.spikes{2};
spikes_adj=spikes;
n_block=0;
for i=1:(length(datarun.triggers)-1)
    actual_t_start=datarun.triggers(i);
    supposed_t_start=n_block*100/120;
    idx1=spikes > actual_t_start;
    idx2=spikes < datarun.triggers(i+1);
    spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
    n_block=n_block+1;
end
spikes_frame_2=floor(spikes_adj*120);


%% BW STA

STA_BW=zeros(32,32,30);
for j=1:length(spikes_frame_2)
    frame_sp=spikes_frame_2(j);
    if frame_sp>=30 && frame_sp<216000
        STA_BW=STA_BW+fitmovie(:,:,(frame_sp-29):frame_sp);
    end
end
STA_BW=STA_BW-min(STA_BW(:));
STA_BW=STA_BW/max(abs(STA_BW(:)));
for i=1:30 
    imagesc(STA_BW(:,:,i)); 
    colormap gray; 
    axis image; 
    pause(0.1); 
end

%% COLOR STA

STA_color=zeros(32,32,3,30);
for j=1:length(spikes_frame)
    frame_sp=spikes_frame(j);
    if frame_sp>=30 && frame_sp<216000
        STA_color=STA_color+fitmovie_color(:,:,:,(frame_sp-29):frame_sp);
    end
end
STA_color=STA_color-min(STA_color(:));
STA_color=STA_color/max(STA_color(:));
for i=1:30
    imagesc(STA_color(:,:,:,i));  
    axis image; 
    pause(0.1); 
end