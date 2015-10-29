function newc=display_sta(index, datarun, hSTA, hZoomSTA, rad)

persistent newCones
global acceptedCones

if ~isempty(newCones)
    acceptedCones = [acceptedCones; newCones];
end
newCones = [];

rad = str2num(rad);

datarun = load_sta(datarun,'load_sta',datarun.cell_ids(index),'keep_java_sta',true);

sta = double(datarun.stas.stas{index}(:,:,:,6));
tmp = abs(sum(sta,3));
sta = sta/max(abs(sta(:)))/2+0.5;
sta = sum(sta,3);
subplot(hSTA)
colormap gray
imagesc(sta)

tmp = tmp/max(tmp(:));

[a,b]=find(tmp>0.7);

subplot(hZoomSTA)
% colormap jet
imagesc(sta)
bord = 10;
axis([min(b)-bord max(b)+bord min(a)-bord max(a)+bord ])

sta = double(datarun.stas.stas{index}(:,:,:,6));

% blue cones have only blue channel; ML cones have all three. Compare
% channels for classifying cones in RGB

if size(sta,3)==3 % RGB
%     rgcones = squeeze(sum(sta(:,:,1:2),3));
%     bcones = squeeze(sta(:,:,3));
    newSTA = mean(sta,3);
else % BW
    newSTA = squeeze(sta(:,:,1));
end

if max(newSTA(:))<abs(min(newSTA(:)))
    newSTA = -newSTA;
end

% estimate SNR
thresh = robust_mean(newSTA(:))+4*robust_std(newSTA(:));

if ~isempty(find(newSTA(:)>thresh, 1))
    
    myMax = max(newSTA(:));
    myInd = 1;
    tmp_sta = newSTA;
    subplot(hZoomSTA)
    hold on    
    while myMax>thresh
        [row, col]=find(tmp_sta==myMax,1);
        tmp_sta(row-rad:row+rad,col-rad:col+rad,:)=0;
        newCones = [newCones; col, row];
        myInd = myInd+1;
        myMax = max(tmp_sta(:));
        plot(col, row,'x', 'markersize',15, 'color','r')
    end    
end
 
        
newc=newCones;




