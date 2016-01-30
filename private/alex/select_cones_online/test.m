rf_size=[9 9];
cone_templ=zeros(9,9,25);
exactCorrection=zeros(25,2);
cnt=1;
for i=-2:2
    for j=-2:2
        bw_kern = make_gaussian('dim',2,'x_size',rf_size(2),'y_size',rf_size(1),'normalize','sum',...
            'center_radius',0.9,'center',[5+i/3,5+j/3], 'effective_radius',2);
        cone_templ(:,:,cnt) = full(bw_kern);
        exactCorrection(cnt,:)=[i/3,j/3];
        cnt=cnt+1;
    end
end

template_contours = cell(25,2);
template_blocks = cell(25,1);
figure
for i=1:25
    b = cone_templ(:,:,i);
    b(b>0.05) = 1;
    b(b<1) = 0;
    dd = imresize(b,5,'method', 'nearest');
    [r, c] = find(dd,1);
    
    [k,m] = find(cone_templ(:,:,i)>0.1);
    
    template_blocks{i} = [k m];
    contour = bwtraceboundary(dd,[r c],'W',8,Inf,'counterclockwise');
    contour= round(contour/5)+0.5;
    subplot(5,5,i)
    imagesc( cone_templ(:,:,i))
    hold on
    plot(contour(:,2), contour(:,1), 'r','linewidth', 2)
    template_contours{i, 1} = contour(:,2);
    template_contours{i, 2} = contour(:,1);

end



figure
colormap gray
imagesc(tmp_sta)


a = tmp_sta;
d_sta = tmp_sta;
figure(22)
colormap gray
imagesc(tmp_sta)
hold on

for kkk=1:5
    myMax=max(d_sta(:));
    [row, col]=find(tmp_sta==myMax,1);
    myCone = tmp_sta(row-2:row+2, col-2:col+2);
    maxCor=zeros(1,25);
    for j=1:25
        A = myCone .* cone_templ(:,:,j);
        maxCor(j) = sum(A(:));
    end
    [~,t]=max(maxCor);
    
    figure(22)
    tp(kkk) = plot(template_contours{t, 1}+col-3, template_contours{t, 2}+row-3, 'r','linewidth', 2)
    
    tk = [template_blocks{t}(:,1)+col-3 template_blocks{t}(:,2)+row-3];
    lininds = zeros(1,length(tk));
    for j=1:length(tk)
        lininds(j) = sub2ind(size(d_sta),tk(j,2),tk(j,1));
    end
    d_sta(lininds) = 0;
end
% figure(23)
% imagesc(d_sta)
% colormap gray

tp(kkk)








tmp = zeros(coarserun.stimulus.field_height, coarserun.stimulus.field_width, length(coarserun.cell_ids));
[columnsInImage rowsInImage] = meshgrid(1:coarserun.stimulus.field_width, 1:coarserun.stimulus.field_height);
for i=1:length(coarserun.cell_ids)
    if ~any(coarserun.stas.fits{i}.sd<0.2)
        sd = mean(coarserun.stas.fits{i}.sd)*1.5;
        m = coarserun.stas.fits{i}.mean;
        
        circlePixels = (rowsInImage - m(2)).^2 + (columnsInImage - m(1)).^2 <= sd.^2;
        
        tmp(:,:,i) = circlePixels;
    end
end


companions = cell(1,length(coarserun.cell_ids));
for i=1:length(coarserun.cell_ids)
    [a,b] = find(tmp(:,:,i));
    allcells = [];
    for j=1:length(a)
        allcells = [allcells; find(tmp(a(j),b(j),:))];
    end
    allcells = unique(allcells);
    allcells(allcells==i) = [];
    companions{i} = allcells;
%     imagesc(tmp(:,:,i))
end

a=companions{datInd};

sta1=datarun.stas.stas{a(1)};
sta1=sta1(:,:,:,params.frame); % frame

sta2=datarun.stas.stas{a(2)};
sta2=sta2(:,:,:,params.frame); % frame

figure
colormap gray
imagesc(sta)


tt = zeros(400,400,3);
tt(:,:,1) = tmp_sta/max(tmp_sta(:));
tt(:,:,2) = sta1/min(sta1(:));
figure
imagesc(tt)


tmp = zeros(400, 400, length(datarun.cell_ids));
for i=1:length(coarserun.cell_ids)
    sta = datarun.stas.stas{i}(:,:,params.frame);
    if abs(min(sta(:)))>max(sta(:)) % OFF cell
        sta = -sta;
    end
    t = robust_std(sta(:))*5;      
    if any(sta(:)>t)
        sta = sta/max(sta(:));
        if sum(sta>0.3)<200
            sta(sta<0.1) = 0;
            tmp(:,:,i) = sta;
            i
        end
    end
end

figure
colormap gray
imagesc(tmp(:,:,datInd))

figure
colormap gray
imagesc(sum(tmp(:,:,[datInd 99]),3))




companions = cell(1,length(datarun.cell_ids));
for i=1:length(datarun.cell_ids)
    [a,b] = find(tmp(:,:,i));
    allcells = [];
    for j=1:length(a)
        allcells = [allcells; find(tmp(a(j),b(j),:))];
    end
    allcells = unique(allcells);
    allcells(allcells==i) = [];
    companions{i} = allcells;
%     imagesc(tmp(:,:,i))
end


tt = zeros(400);
for i=1:400
    for j=1:400
        my_cells = find(tmp(i,j,:)>0.8);
        if ~isempty(my_cells)
            tt(i,j) = 1;
        end
    end
end

figure
colormap gray
imagesc(tt)

figure
colormap gray
imagesc(tmp(:,:,58))



