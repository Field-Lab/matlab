function res = fit_area(datarun, cellID)


sta = squeeze(datarun.stas.stas{cellID});
sta_length = size(sta);
sta_length = sta_length(end);
len = min(sta_length,20)-1;
if size(sta,4)==1
    tc = datarun.vision.timecourses(cellID).g(end-len:end);
    sta = double( sta(:,:,end-len:end));
    sta = sta.*reshape(repmat(tc', size(sta,1)*size(sta,2),1), size(sta,1), size(sta,2),len+1 );
    sta = sum(sta,3);
else
    tcr = datarun.vision.timecourses(cellID).r(end-len:end);
    tcg = datarun.vision.timecourses(cellID).g(end-len:end);
    tcb = datarun.vision.timecourses(cellID).b(end-len:end);
    if ~isempty(tcr) || ~isempty(tcg) || ~isempty(tcb)
        star = double(squeeze( sta(:,:,1,end-len:end)));
        stag = double(squeeze( sta(:,:,2,end-len:end)));
        stab = double(squeeze( sta(:,:,3,end-len:end)));
        star = star.*reshape(repmat(tcr', size(star,1)*size(star,2),1), size(star,1), size(star,2), size(star,3) );
        stag = stag.*reshape(repmat(tcg', size(stag,1)*size(stag,2),1), size(stag,1), size(stag,2), size(stag,3) );
        stab = stab.*reshape(repmat(tcb', size(stab,1)*size(stab,2),1), size(stab,1), size(stab,2), size(stab,3) );
        
        sta = sum(star+stag+stab,3);
    else
        res = 0;
        return;
    end
end

sta = sta/(max(sta(:))*2)+0.5;
sta = imresize(sta, [datarun.stimulus.stixel_height*datarun.stimulus.field_height...
    datarun.stimulus.stixel_width*datarun.stimulus.field_width]);

sta = uint8(round(sta*256));
t = sta([1:30 end-29:end],[1:30 end-29:end])*1.2;
tmp = double(max(t(:)))/255;
bw = im2bw(sta, tmp);
%             bw = im2bw(sta, graythresh(sta));
bw2 = imfill(bw,'holes');

%             figure
%             imshow(bw2)

L = bwlabel(bw2);
tmp = [];
for k = 1:length(unique(L(:)))
    tmp(k) = nnz(L==k);
end
[~, myInd] = max(tmp);
L(L~=myInd) = 0;
L(L==myInd) = 1;
[r,c]=find(L==1,1);
% contour = bwtraceboundary(L,[r c],'W',8,Inf,'counterclockwise');
res = nnz(L);