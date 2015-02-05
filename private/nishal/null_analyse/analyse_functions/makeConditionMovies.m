function makeConditionMovies(cond_str,condMovies,datarun,CellID,CondsToWrite)
for icond=1:2
writerObj = VideoWriter(sprintf('/Volumes/Analysis/nishal/Presentations/Figures for EJ/Null_figures_for_EJ/Movie_cell_%d_cond_%d.avi',CellID,icond));
writerObj.FrameRate=120;

scale=10;
mov=permute(condMovies{CondsToWrite(icond)},[2,3,1]);


mov4D=ones(size(mov,1)*scale,(size(mov,2))*scale,3,size(mov,3)/2 );
h=figure;
X = plot_rf_fit_nps(datarun,CellID,'fill_color',[0,1,0],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[0,0,0]);
centerX = floor(mean(X(1,:)));
centerY = floor(mean(X(2,:)));
X(1,:)=(X(1,:)-centerX)*2+centerX;
X(2,:) = (X(2,:)-centerY)*2+centerY;
X=X*scale;
centerX = floor(mean(X(1,:)));
centerY = floor(mean(X(2,:)));



STA_outline=zeros(size(mov,1)*scale,size(mov,2)*scale);
for iL=1:size(X,2)
STA_outline(end-floor(X(2,iL)),floor(X(1,iL)))=1;    
end

STA_outline=logical(STA_outline);

rows=1:size(STA_outline,1);
cols=1:size(STA_outline,2);
PixelsSelect = zeros(size(STA_outline,1),size(STA_outline,2));
rowsSelect = abs(rows-centerX)<60;
colsSelect = abs(cols-centerY)<60;
PixelsSelect(flip(colsSelect) ,rowsSelect)=1;
PixelsSelect=logical(PixelsSelect);

for itime=1:size(mov,3)/2

xx = imresize(mov(:,:,itime+130),scale,'nearest')'+0.5;
xx(STA_outline)=0;
mov4D(:,1:size(mov,2)*scale,3,itime)=xx;
mov4D(:,1:size(mov,2)*scale,2,itime)=xx;
xx(STA_outline)=1;
mov4D(:,1:size(mov,2)*scale,1,itime)=xx;
end





open(writerObj);
mov4D=mov4D(flip(colsSelect),rowsSelect,:,:);
writeVideo(writerObj,mov4D);

close(writerObj);

%%
end


%% Write STA
h=figure('Color','w');
datarun = load_sta(datarun);
idx=find(CellID==datarun.cell_ids);
sta=datarun.stas.stas{idx};
frame=28;
sta_toUse=zeros(size(sta,1)*scale,size(sta,2)*scale,3);
xx = imresize(sum(sta(:,:,:,frame),3)',scale,'nearest')'+0.5;
xx(STA_outline)=0;
sta_toUse(:,:,3)=xx;
xx = imresize(sum(sta(:,:,:,frame),3)',scale,'nearest')'+0.5;
xx(STA_outline)=0;
sta_toUse(:,:,2)=xx;

xx = imresize(sum(sta(:,:,:,frame),3)',scale,'nearest')'+0.5;
xx(STA_outline)=1;
sta_toUse(:,:,1)=xx;


imagesc(sta_toUse(flip(colsSelect),rowsSelect,:));
caxis([min(sta(:)),max(sta(:))]);
 set(gca,'XColor','w')
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YColor','w');
s=hgexport('readstyle','STA');
hgexport(h,sprintf('/Volumes/Analysis/nishal/Presentations/Figures for EJ/Null_figures_for_EJ/STA_CellID_%d.eps',CellID),s);
  

    
    
    
    
    













end