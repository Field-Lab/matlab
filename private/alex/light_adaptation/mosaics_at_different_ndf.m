%% mosaics at different light levels

%% 2015-03-09-2
% NDF5
% data005 BW-20-10-0.48-11111-16x16 1800s
% data006 NSEM repeats
% data007 BW-20-10-0.48-11111 repeats

% NDF4
% data008 BW-10-8-0.48-11111 1800s
% data009 NSEM repeats
% data010 BW-10-8-0.48-11111 repeats

% NDF3
% data011 BW-16-8-0.48-11111 1800s
% data012 NSEM repeats
% data013 BW-8-8-0.48-11111 repeats
% data014 BW-8-8-0.48-11111 1386s

% NDF2
% data015 BW-16?6-0.48-11111 900s
% data016 NSEM repeats
% data017 BW-8-6-0.48-11111 repeats
% data018 BW-8-6-0.48-11111 1200s

% NDF1
% data019 RGB-16-6-0.48-11111 900s
% data020 NSEM repeats
% data021 BW-8-6-0.48-11111 repeats
% data022 BW-8-6-0.48-11111 1200s

% NDF0
% data023 RGB-16-4-0.48-11111 900s
% data024 NSEM repeats
% data025 BW-4-2-0.48-11111 repeats
% data026 BW-4-2-0.48-11111 1800s
% data027 RGB-4-2-0.48-11111 1800s


path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/';

staruns=['data005'; 'data008'; 'data011'; 'data014'; 'data015'; 'data018'; ...
    'data019'; 'data022'; 'data023'; 'data026'; 'data027'];

ndfs = [5 4 3 3 2 2 1 1 0 0 0];
scales = [20 10 16 8 16 8 16 8 16 4 4];

% % file path to save pictures
% filepath=['/Users/alexth/Desktop/Light_adaptation/NSEM_WN/',date,'/',concatname,'/'];
% if ~exist(filepath,'dir')
%     mkdir(filepath);
% end

for i=1:size(staruns,1)
    datarun = load_data(fullfile(path2data, staruns(i,:), staruns(i,:)));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = get_sta_fits_from_vision(datarun); 
%     datarun = load_sta(datarun);
%     datarun = set_polarities(datarun);

    figure
    set(gcf,'Name',['NDF ', int2str(ndfs(i))])
    cnt = 1;
    for j=[1:5 8 9 11 12]
        subplot(3,3,cnt)
        plot_rf_summaries(datarun, {j}, 'clear', false, 'plot_fits', true, 'fit_color', 'k')
        title(datarun.cell_types{j}.name)
        cnt=cnt+1;
    end
end

% compare mosaics with different STA run
% NDF0
colors = 'kbr';
cnt_colors=1;
for i=9:11
    datarun = load_data(fullfile(path2data, staruns(i,:), staruns(i,:)));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = get_sta_fits_from_vision(datarun); 

    cnt = 1;
    for j=[1:5 8 9 11 12]
        subplot(3,3,cnt)
        plot_rf_summaries(datarun, {j}, 'clear', false, 'scale', scales(i), 'plot_fits', true, 'fit_color', colors(cnt_colors))
        title(datarun.cell_types{j}.name)
        cnt=cnt+1;
    end
    cnt_colors = cnt_colors+1;
end


for_use = ['data005'; 'data008'; 'data014'; 'data018'; 'data022'; 'data026'];

scales = [20 10 8 8 8 4];
colors = 'kbrmcy';
cnt_colors=1;
for i=2:size(for_use,1)-1
    datarun = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = get_sta_fits_from_vision(datarun); 
    
    figure

    cnt = 1;
    for j=1:4%[1:5 8 9 11 12]
        
        subplot(2,2,cnt)
        plot_rf_summaries(datarun, {j}, 'clear', false, 'scale', scales(i), 'plot_fits', true, 'fit_color', colors(cnt_colors))
        title(datarun.cell_types{j}.name)
        cnt=cnt+1;
    end
    cnt_colors = cnt_colors+1;
    drawnow
    
    datarun = load_data(fullfile(path2data, for_use(i+1,:), for_use(i+1,:)));
    datarun = load_params(datarun,'verbose',1);
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = get_sta_fits_from_vision(datarun); 
    
    cnt = 1;
    for j=1:4%[1:5 8 9 11 12]
        
        subplot(2,2,cnt)
        plot_rf_summaries(datarun, {j}, 'clear', false, 'scale', scales(i+1), 'plot_fits', true, 'fit_color', colors(cnt_colors))
        title(datarun.cell_types{j}.name)
        cnt=cnt+1;
    end
    cnt_colors = cnt_colors-1;
    drawnow
    
    
end


%% cell by cell basis - mosaics from combined RFs

refresh_times = [10 8 8 6 6 2];
datarun = cell(6,1);
for i=1:6
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_sta(datarun{i},'load_sta',[]);
    datarun{i} = get_sta_fits_from_vision(datarun{i});
    datarun{i} = load_sta(datarun{i},'load_sta', 'all');
end

col = 'rcbmkg';
t = zeros(6,100,4);
time_to_peak =zeros(6,100,4);
half_width = zeros(6,100,4);
for cell_type = 1
    cellIDs = get_cell_indices(datarun{1}, {cell_type});
    
    figure
    set(gcf,'Name',datarun{1}.cell_types{cell_type}.name)
    
    for i=1:6 % light levels
        
        subplot(2,3,i)
        
        for j=1:length(cellIDs)
            
            sta = squeeze(datarun{i}.stas.stas{cellIDs(j)});
            tc = datarun{i}.vision.timecourses(cellIDs(j)).g(11:end);
            for_tc = imresize(sta, scales(i));
            sta = sta(:,:,11:end);            
            sta = sta.*reshape(repmat(tc', size(sta,1)^2,1), size(sta,1), size(sta,2),20 );
            sta = sum(sta,3);
            sta = sta/(max(sta(:))*2)+0.5;
            sta = imresize(sta, scales(i));
            sta = uint8(round(sta*256));
            
            r = sta([1:30 290:end],[1:30 290:end])*1.2;
            tmp = double(max(r(:)))/255;
            bw = im2bw(sta, tmp);           
%             bw = im2bw(sta, graythresh(sta));

            bw2 = imfill(bw,'holes');
            L = bwlabel(bw2);
            
            tmp = [];
            for k = 1:length(unique(L(:)))
                tmp(k) = nnz(L==k);
            end
            if tmp(1)~=0
                [~, myInd] = max(tmp);
                L(L~=myInd) = 0;
                L(L==myInd) = 1;
                [r,c]=find(L==1,1);
                contour = bwtraceboundary(L,[r c],'W',8,Inf,'counterclockwise');
                if isempty(find(r<30,1)) && isempty(find(r>290,1)) && isempty(find(c<30,1)) && isempty(find(c>290,1))
                    
                    my_tc = 0;
                    
                    if cell_type<3
                        if size(contour,1)<200
                            hold on;
%                             plot(contour(:,2),contour(:,1),col(mod(j,length(col))+1),'LineWidth',1);
                            t(i,j, cell_type) = nnz(L);
                            
                            for k = 1:size(contour,1)
                                my_tc = my_tc+squeeze(for_tc(contour(k,1),contour(k,2),:));
                            end                            
                        time_axis = -linspace(0, refresh_times(i)*(1000/120)*29,30);
                         hold on
                        plot(time_axis, my_tc(end:-1:1)/max(abs(my_tc)));
                        z = find(abs(my_tc(end:-1:1))/max(abs(my_tc))==1, 1);
                        time_to_peak(i,j, cell_type) = time_axis(z);
                            
                        end
                    elseif size(contour,1)<150
                        hold on;
%                         plot(contour(:,2),contour(:,1),col(mod(j,length(col))+1),'LineWidth',1);
                        t(i,j, cell_type) = nnz(L);
                        
                        for k = 1:size(contour,1)
                            my_tc = my_tc+squeeze(for_tc(contour(k,1),contour(k,2),:));
                        end
                        time_axis = -linspace(0, refresh_times(i)*(1000/120)*29,30);
                         hold on
                        plot(time_axis, my_tc(end:-1:1)/max(abs(my_tc)));
                        z = find(abs(my_tc(end:-1:1))/max(abs(my_tc))==1, 1);
                        time_to_peak(i,j, cell_type) = time_axis(z);
                        
                    end
                    if length(my_tc)>1
                        new_tc = time_axis(end):1/3:0;
                        a = interp1(time_axis(end:-1:1),my_tc/max(abs(my_tc)), new_tc);                        
                        a = abs(a(end:-1:1));
                        [~, ind]=max(a);
                        beg = find(a>0.5,1);
                        ending = find(a(beg:end)<0.5,1)+beg-1;
                        half_width(i,j,cell_type) = abs(new_tc(end - ending)-new_tc(end - beg));
                    end
                end
            end
     
        end
        axis([-1500 0 -1 1])
        axis off
%         axis ij
%         set(gca, 'Dataaspectratio', [1 1 1])
%         axis([1 320 1 320]);
%         axis off
        
    end
    
end


pixel_size = 5.5; % mkm

mean_area = zeros(4,6);
mean_std = zeros(4,6);
mean_peak = zeros(4,6);
mean_peak_std = zeros(4,6);
for i=1:4
    for j=1:6
        tmp = t(j,:,i);
        r = find(tmp>=100&tmp<8000);
        tmp = tmp(r);
        mean_area(i,j) = round(robust_mean(tmp)*(pixel_size^2));
        mean_std(i,j) = round(robust_std(tmp)*(pixel_size^2)/sqrt(length(tmp)));
        
        tmp = half_width(j,r,i);

        mean_peak(i,j) = robust_mean(tmp);
        mean_peak_std(i,j) = robust_std(tmp)/sqrt(length(tmp));        
    end
end

figure
errorbar(mean_area', mean_std')
legend('ON P', 'OFF P', 'ON m', 'OFF m')
set(gca,'xtick', 1:6, 'xticklabel', {'5','4', '3', '2', '1', '0'})
xlabel('NDF')
ylabel('area, mkm^2')
title('Absolute RF area')


figure
errorbar(mean_peak', mean_peak_std')
legend('ON P', 'OFF P', 'ON m', 'OFF m', 'location', 'best')
set(gca,'xtick', 1:6, 'xticklabel', {'5','4', '3', '2', '1', '0'})
xlabel('NDF')
ylabel('half-height width, ms')
title('Time course')


figure
plot((mean_area./repmat(mean_area(:,end),1,6))')
legend('ON P', 'OFF P', 'ON m', 'OFF m')
set(gca,'xtick', 1:6, 'xticklabel', {'5','4', '3', '2', '1', '0'})
xlabel('NDF')
ylabel('area, relative to NDF0')
title('Relative RF area')


%% get example STA and time course


refresh_times = [10 8 8 6 6 2];
subpl = [1 2 3 5 6 7];
for cell_type = 1:4
    cellIDs = get_cell_indices(datarun{1}, {cell_type});
    
    for j=3%length(cellIDs)
        
        figure
        set(gcf,'Name',[datarun{1}.cell_types{cell_type}.name, '  cell ID ', int2str(cellIDs(j))])
        set(gcf, 'position', [-1522         381        1092         724])
        
        for i=1:6 % light levels
            
            subplot(2,4,subpl(i))
            sta = squeeze(datarun{i}.stas.stas{cellIDs(j)});
            tc = datarun{i}.vision.timecourses(cellIDs(j)).g(11:end);
            sta = sta(:,:,11:end);            
            sta = sta.*reshape(repmat(tc', size(sta,1)^2,1), size(sta,1), size(sta,2),20 );
            sta = sum(sta,3);
            sta = sta/(max(sta(:))*2)+0.5;
            sta = imresize(sta, scales(i));
            
            colormap gray
            imshow(1-sta)
            hold on
            
            sta = uint8(round(sta*256));
            t = sta([1:30 290:end],[1:30 290:end])*1.2;
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
            contour = bwtraceboundary(L,[r c],'W',8,Inf,'counterclockwise');
            if isempty(find(r<30,1)) && isempty(find(r>290,1)) && isempty(find(c<30,1)) && isempty(find(c>290,1))
                
                if cell_type<3
                    if size(contour,1)<200
                        hold on;
                        plot(contour(:,2),contour(:,1),'r','LineWidth',1);
                        t(i,j, cell_type) = nnz(L);
                        subplot(2,4,8)
                        hold on;
                        plot(contour(:,2),contour(:,1),'linewidth',2)
                    end
                elseif size(contour,1)<150
                    hold on;
                    plot(contour(:,2),contour(:,1),'r','LineWidth',1);
                    t(i,j, cell_type) = nnz(L);
                    subplot(2,4,8)
                    hold on;
                    plot(contour(:,2),contour(:,1), 'linewidth',2)
                end
            end
            
            % get new time course
            sta = squeeze(datarun{i}.stas.stas{cellIDs(j)});
            sta = imresize(sta, scales(i)).*repmat(L,1,1,30);
            tc = sum(sum(sta));
            tc = tc/sum(L(:));
            time_axis = -linspace(0, refresh_times(i)*(1000/120)*19,20);
            subplot(2,4,4)
            hold on
            plot(time_axis, squeeze(tc(:,:,end:-1:11)), 'linewidth',2)            
        end
        
        
        subplot(2,4,4)
        axis tight
        set(gca, 'xtick', [-1200:300:0], 'xticklabel', cellfun(@int2str,{([-1200:300:0])'}, 'UniformOutput', 0))
        legend('5', '4', '3', '2', '1', '0',  'location', 'best')
        line([-1600 0], [0 0], 'color','k')
        set(gca, 'XLim', [-1000, 0])
        
        subplot(2,4,8)
        axis tight
        axis ij
        for k= subpl
            subplot(2,4,k)
            set(gca, 'Dataaspectratio', [1 1 1])
        end
        
    end
    
end
axis off






%% cell by cell basis - mosaics from combined RFs

for_use = ['data006';'data007'; 'data011'];

datarun = cell(3,1);
for i=1:3
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_sta(datarun{i},'load_sta',[]);
    datarun{i} = get_sta_fits_from_vision(datarun{i});
    datarun{i} = load_sta(datarun{i},'load_sta', 'all');
end

scales = [16 16 10];
col = 'rcbmkg';
t = zeros(3,100,4);
for cell_type = 1:4
    cellIDs = get_cell_indices(datarun{1}, {cell_type});
    
    figure
    set(gcf,'Name',datarun{1}.cell_types{cell_type}.name)
    
    for i=1:3 % light levels
        
        subplot(1,3,i)
        
        for j=1:length(cellIDs)
            
            sta = squeeze(datarun{i}.stas.stas{cellIDs(j)});
            tc = datarun{i}.vision.timecourses(cellIDs(j)).g(11:end);
            sta = sta(:,:,11:end);
            sta = sta.*reshape(repmat(tc', size(sta,1)^2,1), size(sta,1), size(sta,2),20 );
            sta = sum(sta,3);
            sta = sta/(max(sta(:))*2)+0.5;
            sta = imresize(sta, scales(i));
            sta = uint8(round(sta*256));
            
            t = sta([1:30 290:end],[1:30 290:end])*1.2;
            tmp = double(max(t(:)))/255;
            bw = im2bw(sta, tmp);           
%             bw = im2bw(sta, graythresh(sta));

            bw2 = imfill(bw,'holes');
            L = bwlabel(bw2);
            
            tmp = [];
            for k = 1:length(unique(L(:)))
                tmp(k) = nnz(L==k);
            end
            [~, myInd] = max(tmp);
            L(L~=myInd) = 0;
            L(L==myInd) = 1;
            [r,c]=find(L==1,1);
          
            contour = bwtraceboundary(L,[r c],'W',8,Inf,'counterclockwise');
            if isempty(find(r<30,1)) && isempty(find(r>290,1)) && isempty(find(c<30,1)) && isempty(find(c>290,1))
                
                if cell_type<3
                    if size(contour,1)<200
                        hold on;
                        plot(contour(:,2),contour(:,1),col(mod(j,length(col))+1),'LineWidth',2);
                        t(i,j, cell_type) = nnz(L);
                    end
                elseif size(contour,1)<150
                    hold on;
                    plot(contour(:,2),contour(:,1),col(mod(j,length(col))+1),'LineWidth',2);
                    t(i,j, cell_type) = nnz(L);
                end
            end
     
        end
        set(gca, 'Dataaspectratio', [1 1 1])
        axis([1 320 1 320]);
        
    end
    
end




pixel_size = 5.5; % mkm

mean_area = zeros(4,3);
mean_std = zeros(4,3);
for i=1:4
    for j=1:3
        tmp = t(j,:,i);
        tmp = tmp(tmp>=100&tmp<8000);
        mean_area(i,j) = round(robust_mean(tmp)*(pixel_size^2));
        mean_std(i,j) = round(robust_std(tmp)*(pixel_size^2)/sqrt(length(tmp)));
    end
end

figure
errorbar(mean_area', mean_std')
legend('ON P', 'OFF P', 'ON m', 'OFF m')
set(gca,'xtick', 1:3, 'xticklabel', {'4', '3', '2'})
xlabel('NDF')
ylabel('area, mkm^2')
title('Absolute RF area')

figure
plot((mean_area./repmat(mean_area(:,end),1,3))')
legend('ON P', 'OFF P', 'ON m', 'OFF m')
set(gca,'xtick', 1:3, 'xticklabel', {'4', '3', '2'})
xlabel('NDF')
ylabel('area, relative to NDF0')
title('Relative RF area')



%% All Mosaics at each light level

datarun = cell(6,1);
for i=1:6
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_sta(datarun{i},'load_sta',[]);
    datarun{i} = get_sta_fits_from_vision(datarun{i});
    datarun{i} = load_sta(datarun{i},'load_sta', 'all');
end

for    i = 2:6;% NDF0
    
    figure
    col = 'k'%'rcbmkg';
    t = zeros(6,100,4);
    for cell_type = 1:4
        cellIDs = get_cell_indices(datarun{1}, {cell_type});
        
        %     figure
        %     set(gcf,'Name',datarun{1}.cell_types{cell_type}.name)
        
        
        subplot(2,2,cell_type)
        
        for j=1:length(cellIDs)
            
            sta = squeeze(datarun{i}.stas.stas{cellIDs(j)});
            tc = datarun{i}.vision.timecourses(cellIDs(j)).g(11:end);
            sta = sta(:,:,11:end);
            for_tc = imresize(sta, scales(i));
            sta = sta.*reshape(repmat(tc', size(sta,1)^2,1), size(sta,1), size(sta,2),20 );
            sta = sum(sta,3);
            sta = sta/(max(sta(:))*2)+0.5;
            sta = imresize(sta, scales(i));
            sta = uint8(round(sta*256));
            
                        
            t = sta([1:30 290:end],[1:30 290:end])*1.2;
            tmp = double(max(t(:)))/255;
            bw = im2bw(sta, tmp);           
%             bw = im2bw(sta, graythresh(sta));

            bw2 = imfill(bw,'holes');
            L = bwlabel(bw2);
            
            tmp = [];
            for k = 1:length(unique(L(:)))
                tmp(k) = nnz(L==k);
            end
            if tmp(1)~=0
                [~, myInd] = max(tmp);
                L(L~=myInd) = 0;
                L(L==myInd) = 1;
                [r,c]=find(L==1,1);
                contour = bwtraceboundary(L,[r c],'W',8,Inf,'counterclockwise');
                if isempty(find(r<30,1)) && isempty(find(r>290,1)) && isempty(find(c<30,1)) && isempty(find(c>290,1))
                    
                    my_tc = 0;
                    if cell_type<3
                        if size(contour,1)<200
                            hold on;
%                             plot(contour(:,2),contour(:,1),col(mod(j,length(col))+1),'LineWidth',2);
                            t(i,j, cell_type) = nnz(L);
                            for k = 1:size(contour,1)
                                my_tc = my_tc+squeeze(for_tc(contour(k,1),contour(k,2),:));
                            end
                            plot(my_tc/max(abs(my_tc)));
                        end
                    elseif size(contour,1)<150
                        hold on;
%                         plot(contour(:,2),contour(:,1),col(mod(j,length(col))+1),'LineWidth',2);
                        t(i,j, cell_type) = nnz(L);
                        for k = 1:size(contour,1)
                            my_tc = my_tc+squeeze(for_tc(contour(k,1),contour(k,2),:));
                        end
                        plot(my_tc/max(abs(my_tc)));
                    end

                end
            end
            
        end
%         set(gca, 'Dataaspectratio', [1 1 1])
%         axis([1 320 1 320]);
        axis off
        
    end
end

