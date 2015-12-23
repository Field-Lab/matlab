
temp = load('~/Dropbox/Lab/Development/matlab-standard/private/freddy/512elecpositions.mat'); % Find a more general location for this or call a different text file.
positions = temp.positions;

load('/Volumes/Analysis/nishal/data/all_1min.mat');
%% Raw movie
f = figure; set(f,'Position',[100 360 1000 550]);
set(f,'Color','white');
for itime=1000:1150
    
        scatter(positions(:,1),positions(:,2),350,data(:,itime),'filled'); 
        axis off; axis image; colorbar;
        caxis([-300 300]); 
        
        hold on; 
       % scatter(positions(stimChan,1),positions(stimChan,2),350,'black');

    
scatter(positions(:,1),positions(:,2),350,data(:,itime),'black');    
pause(0.001);
end

%% 
near_elecs=cell(512,1);
for i=1:512
near_elecs{i} = find_near_elecs(positions, i);
end

figure;
ielec=512;
scatter(positions(:,1),positions(:,2),350,data(:,itime),'blue');    
hold on 
scatter(positions(near_elecs{ielec},1),positions(near_elecs{ielec},2),350,data(near_elecs{ielec},itime),'black');    

%% 
data_mean = mean(data,2);
time_win=10;
icnt=0;
data_rms=zeros(512,1);
for itime=1000:1500
    icnt=icnt+1;
    
    for ielec=1:512
        da=data(near_elecs{ielec},itime-time_win:itime+time_win)-repmat(data_mean(near_elecs{ielec}),1,(2*time_win+1));
        data_rms(ielec,icnt)=norm(da(:));
    end
end

%%
f = figure; set(f,'Position',[100 360 1000 550]);
set(f,'Color','white');

for iicnt=1:icnt
        scatter(positions(:,1),positions(:,2),350,data_rms(:,iicnt),'filled'); 
        axis off; axis image; colorbar;
        %caxis([-300 300]); 
        
        hold on; 
       % scatter(positions(stimChan,1),positions(stimChan,2),350,'black');

    
scatter(positions(:,1),positions(:,2),350,data_rms(:,iicnt),'black');    
title(sprintf('Time: %0.005f ms',iicnt/20000));
pause(0.001);
end


