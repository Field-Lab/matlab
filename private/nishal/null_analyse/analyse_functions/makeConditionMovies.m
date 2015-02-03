function makeConditionMovies(cond_str,condMovies,datarun,CellID,CondsToWrite)

writerObj = VideoWriter(sprintf('/Volumes/Analysis/nishal/Presentations/Figures for EJ/Null_figures_for_EJ/Movie_cell_%d.avi',CellID));
writerObj.FrameRate=120;

scale=10;
mov=permute(condMovies{CondsToWrite(1)},[2,3,1]);


mov4D=ones(size(mov,1)*scale,(size(mov,2)*2.1)*scale,3,size(mov,3));
h=figure;
X = plot_rf_fit_nps(datarun,CellID,'fill_color',[0,1,0],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[0,0,0]);
X=X*scale;
STA_outline=zeros(size(mov,1)*scale,size(mov,2)*scale);
for iL=1:size(X,2)
STA_outline(end-floor(X(2,iL)),floor(X(1,iL)))=1;    
end

STA_outline=logical(STA_outline);



for itime=1:size(mov,3)

xx = imresize(mov(:,:,itime),scale,'nearest')'+0.5;
xx(STA_outline)=0;
mov4D(:,1:size(mov,2)*scale,3,itime)=xx;
mov4D(:,1:size(mov,2)*scale,2,itime)=xx;
xx(STA_outline)=1;
mov4D(:,1:size(mov,2)*scale,1,itime)=xx;
end


mov=permute(condMovies{CondsToWrite(2)},[2,3,1]);

for itime=1:size(mov,3)
xx = imresize(mov(:,:,itime),scale,'nearest')'+0.5;
xx(STA_outline)=0;
mov4D(:,end-size(mov,2)*scale+1:end,3,itime)=xx;
mov4D(:,end-size(mov,2)*scale+1:end,2,itime)=xx;
xx(STA_outline)=1;
mov4D(:,end-size(mov,2)*scale+1:end,1,itime)=xx;
end




open(writerObj);

writeVideo(writerObj,mov4D);

close(writerObj);

%%


end