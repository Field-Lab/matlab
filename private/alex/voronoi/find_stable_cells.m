
thresh = 0.4;

ww=[]; pixn = []; snr = [];

figure
for i=1:length(datarun.cell_types{4}.cell_ids)
    visionID = datarun.cell_types{4}.cell_ids(i);
    cellInd = find(datarun.cell_ids == visionID);
    
    sta = double(squeeze(datarun.stas.stas{cellInd}));
    sta = sta/min(sta(:));
    [a, b] = find(sta(:,:,4)>0.5);
    tc=0;
    for i=1:length(a)
        tc = tc + squeeze(sta(a(i),b(i),:));
    end
    tc = tc/max(tc);
    tc = repmat(tc, 1, 300);
    tc = repmat(tc, 1, 1, 300);
    tc = shiftdim(tc, 1);
    
    sta1 = sta .* tc;
    sta1 = sum(sta1,3);
    sta1 = sta1/max(sta1(:));
    [a1, b1] = find(sta1>=thresh);
    
    
    sta = double(squeeze(datarun1.stas.stas{cellInd}));
    sta = sta/min(sta(:));
    [a, b] = find(sta(:,:,4)>0.5);
    tc=0;
    for j=1:length(a)
        tc = tc + squeeze(sta(a(j),b(j),:));
    end
    tc = tc/max(tc);
    tc = repmat(tc, 1, 300);
    tc = repmat(tc, 1, 1, 300);
    tc = shiftdim(tc, 1);
    
    sta1 = sta .* tc;
    sta1 = sum(sta1,3);
    sta1 = sta1/max(sta1(:));
    [a2, b2] = find(sta1>=thresh);
    sn = robust_std(sta1(:));
    %     [~,ic]=sort(sta1(:));
%     [a2, b2] = find(sta1>sta1(ic(end - length(a1))));
    
    
%     figure
%     colormap gray
%     imagesc(sta1)
%     hold on
%     plot(b1,a1,'xr')
%     plot(b2,a2,'+b')
    
    
    
    X = [b1, a1];
    Y = [b2, a2];
    
    opt.method='nonrigid';
    [Transform, C]=cpd_register(X,Y, opt);
%     tmp = Y - Transform.Y;
    tmp = pdist2(Y', Transform.Y');
    ww = [ww tmp(2)];
    pixn = [pixn length(a1)];
    snr = [snr sn];
    
end
% 
% figure,cpd_plot_iter(X, Y); title('Before');
% figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');



figure
plot(ww./pixn)
hold on
figure
plot(pixn)
plot(snr)
hold on
plot(ww/1000)


% find cones
mycone = 0;
all_sta = zeros(300,300,length(datarun.cell_types{4}.cell_ids));
for i=1:length(datarun.cell_types{4}.cell_ids)
    visionID = datarun.cell_types{4}.cell_ids(i);
    cellInd = find(datarun.cell_ids == visionID);
    
    sta = double(squeeze(datarun.stas.stas{cellInd}));
    sta = sta/min(sta(:));
    [a, b] = find(sta(:,:,4)>0.5);
    tc=0;
    for j=1:length(a)
        tc = tc + squeeze(sta(a(j),b(j),:));
    end
    tc = tc/max(tc);
    tc = repmat(tc, 1, 300);
    tc = repmat(tc, 1, 1, 300);
    tc = shiftdim(tc, 1);
    
    sta1 = sta .* tc;
    sta1 = sum(sta1,3);
    sta1 = sta1/max(sta1(:));
    all_sta(:,:,i) = sta1;
    [a1, b1] = find(sta1>=thresh);
    
%     figure
%     colormap gray
%     imagesc(sta1)
    [a,b] = find(sta1==1);
    mycone = mycone + sta1(a-10:a+10, b-10:b+10);
end
mycone = mycone/max(mycone(:));
figure
imagesc(mycone)


for i=1:length(datarun.cell_types{4}.cell_ids)
    visionID = datarun.cell_types{4}.cell_ids(i);
    cellInd = find(datarun.cell_ids == visionID);
    
    sta = all_sta(:,:,i);
    sta = sta/max(sta(:));
    sta1=sta;
    mc = 0;
    cnt=1;
    clear a b cw
    tt=robust_std(sta1(:));
    while max(sta(:))>tt*5
        sta1 = sta1/max(sta1(:));
        [a(cnt),b(cnt)] = find(sta1==1);
        mc = mc + sta(a(cnt)-2:a(cnt)+2,b(cnt)-2:b(cnt)+2);
        cw(cnt) = sta(a(cnt),b(cnt));
        sta(a(cnt)-2:a(cnt)+2,b(cnt)-2:b(cnt)+2) = 0;
        sta1 = sta;
        cnt = cnt+1;
    end
%     figure
%     sta = all_sta(:,:,i);
%     colormap gray
%     imagesc(sta)
%     hold on
%     plot(b,a,'rx')
    sta_cones{i} = [b;a; cw];
end




% find cones STA1
mycone = 0;
all_sta = zeros(300,300,length(datarun.cell_types{4}.cell_ids));
for i=1:length(datarun.cell_types{4}.cell_ids)
   
    visionID = datarun1.cell_types{4}.cell_ids(i);
    cellInd = find(datarun1.cell_ids == visionID);
    
    sta = double(squeeze(datarun1.stas.stas{cellInd}));
    sta = sta/min(sta(:));
    [a, b] = find(sta(:,:,4)>0.5);
    tc=0;
    for j=1:length(a)
        tc = tc + squeeze(sta(a(j),b(j),:));
    end
    tc = tc/max(tc);
    tc = repmat(tc, 1, 300);
    tc = repmat(tc, 1, 1, 300);
    tc = shiftdim(tc, 1);
    
    sta1 = sta .* tc;
    sta1 = sum(sta1,3);
    sta1 = sta1/max(sta1(:));
    all_sta(:,:,i) = sta1;
    [a1, b1] = find(sta1>=thresh);
    
%     figure
%     colormap gray
%     imagesc(sta1)
    [a,b] = find(sta1==1);
    mycone = mycone + sta1(a-10:a+10, b-10:b+10);
end
mycone = mycone/max(mycone(:));
figure
imagesc(mycone)


% find cones
for i=1:length(datarun.cell_types{4}.cell_ids)
    
    visionID = datarun1.cell_types{4}.cell_ids(i);
    cellInd = find(datarun1.cell_ids == visionID);
   
    
    sta = all_sta(:,:,i);
    sta = sta/max(sta(:));
    sta1=sta;
    mc = 0;
    cnt=1;
    clear a b cw
    tt=robust_std(sta1(:));
    while max(sta(:))>tt*5
        sta1 = sta1/max(sta1(:));
        [a(cnt),b(cnt)] = find(sta1==1);
        mc = mc + sta(a(cnt)-2:a(cnt)+2,b(cnt)-2:b(cnt)+2);
        cw(cnt) = sta(a(cnt),b(cnt));
        sta(a(cnt)-2:a(cnt)+2,b(cnt)-2:b(cnt)+2) = 0;
        sta1 = sta;
        cnt = cnt+1;
    end
%     figure
%     sta = all_sta(:,:,i);
%     colormap gray
%     imagesc(sta)
%     hold on
%     plot(b,a,'rx')
    sta_cones1{i} = [b;a; cw];
end



cnt=1;
ww=[];
for i=29:44%length(datarun.cell_types{4}.cell_ids)
    
   a = sta_cones{i}(3,:);
   b = sta_cones1{i}(3,:);
   if length(a)>length(b)
       sta_cones{i} = sta_cones{i}(:,1:length(b));
   else
       sta_cones1{i} = sta_cones1{i}(:,1:length(a));
   end
    
    X = sta_cones{i}(1:2,:)';
    Y = sta_cones1{i}(1:2,:)';
    
    figure(1)
    subplot(4,4,cnt)
    opt.method='nonrigid';
    [Transform, C]=cpd_register(X,Y, opt);
    tmp = Y - Transform.Y;
%     tmp = pdist2(Y', Transform.Y');
    ww = [ww sum(abs(tmp))];
    drawnow
    pause(0.5)
    
    figure(2)
    subplot(4,4,cnt)
    plot(sta_cones{i}(1,:), sta_cones{i}(2,:), 'r*')
    hold on
    plot(sta_cones1{i}(1,:), sta_cones1{i}(2,:), 'bo')
    title([num2str(ww(end)), ', i=', int2str(i)])
    drawnow
    pause(0.5)
    
    cnt=cnt+1;
end




