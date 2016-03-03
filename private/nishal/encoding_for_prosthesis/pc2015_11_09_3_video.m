

figure('Color','w'); 
plot_rf_fit(datarun,'ON Parasol','fill_color',[0,1,1],'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
axis equal
xlim([15.5,65]);ylim([8,37]);
set(gca,'visible','off');


figure('Color','w'); 
plot_rf_fit(datarun,'OFF Parasol','fill_color',[1,1,0],'fill',false,'edge',true,'labels',false,'edge_color',[0,0,1]);
axis equal
set(gca,'visible','off');
xlim([15.5,65]);ylim([8,37]);

%% 
  
ielec = 81;
files = ls(sprintf('/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/elecRespAuto_n*_p%d.mat',ielec));
files = strsplit(files,'\n');
    
data=cell(length(files),1);
for ifile = 1:length(files)-1
data{ifile}=load(files{ifile});    
end

myVideo = VideoWriter('/Volumes/Lab/Users/bhaishahster/QIF/stimulation2.avi', 'Uncompressed AVI');

myVideo.FrameRate = 15;  % Default 30
%myVideo.Quality = 100;    % Default 75
open(myVideo);

nCell=length(files)-1;
for iAmp=1:39;
    close all
h=figure;
plot_rf_fit(datarun,'ON Parasol','fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
hold on;
plot_rf_fit(datarun,'OFF Parasol','fill',false,'edge',true,'labels',false,'edge_color',[0,0,1]);
axis equal

set(gca,'visible','off');
hold on;
xlim([15.5,65]);ylim([8,37]);
for icell=1:nCell
    color = data{icell}.elecRespAuto.LogisticReg(iAmp);
    cellID = data{icell}.elecRespAuto.neuronInfo.neuronIds;
    plot_rf_fit(datarun,cellID,'fill_color',[0,0,0],'alpha',color,'fill',true,'edge',false,'labels',false);
    axis equal
    set(gca,'visible','off');
    hold on;
    xlim([15.5,65]);ylim([8,37]);
       set(gca,'LooseInset',get(gca,'TightInset'))
       title(sprintf('electrode: %d, current: %f A',ielec,data{icell}.elecRespAuto.stimInfo.listAmps(iAmp)));
end
F = getframe(h);
writeVideo(myVideo, F);

end
close(myVideo);

