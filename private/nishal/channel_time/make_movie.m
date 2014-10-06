
temp = load('~/Dropbox/Lab/Development/matlab-standard/private/freddy/512elecpositions.mat'); % Find a more general location for this or call a different text file.
positions = temp.positions;

load('/Volumes/Analysis/nishal/data/all_1min.mat');
%%
f = figure; set(f,'Position',[100 360 1000 550]);
set(f,'Color','white');
for itime=1:50
    
        scatter(positions(:,1),positions(:,2),350,data(:,itime),'filled'); 
        axis off; axis image; colorbar;
        caxis([-300 300]); 
        
        hold on; 
       % scatter(positions(stimChan,1),positions(stimChan,2),350,'black');

    
scatter(positions(:,1),positions(:,2),350,data(:,itime),'black');    
pause(0.001);
end

