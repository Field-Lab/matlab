/Volumes/Lab/Development/scripts/grind -p -c /Volumes/Lab/Development/vision-xml/current/primate-1cone_ath.xml /Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11 BW-2-6-0.48-22222-300x300-60.35.xml
/Volumes/Lab/Development/scripts/grind -p -c /Volumes/Lab/Development/vision-xml/current/primate-1cone_ath.xml /Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11 BW-2-6-0.48-11111-300x300-60.35.xml

% BW-2-6-0.48-11111-300x300-60.35

clear

local_path = '/Volumes/Analysis/';

datarun = load_data([local_path, '2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11']);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = set_polarities(datarun);
datarun = load_cones_ath(datarun, 'd08-11-norefit_d08-bayes-msf_5.00');
datarun = make_mosaic_struct(datarun);
cell_type = 13;

datarun.cell_types{cell_type}.cell_ids
cellIDs = get_cell_indices(datarun, {cell_type});

datarun = load_sta(datarun, 'load_sta', {cell_type});

sta = squeeze(datarun.stas.stas{112});

figure
colormap gray
imagesc(sta(:,:,4))

figure
set(gcf, 'position', [ -916   484   422   423]);
hold on
plot_rf_summaries(datarun, {cell_type}, 'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'r')
axis([130,210,70,155])
set(gca,'dataaspectratio', [1 1 1])
axis off
axis ij

load('/Volumes/Analysis/2011-12-13-2/cone_data/d08-11-norefit_d08-bayes-msf_5.00/mydll.mat')

figure
set(gcf, 'position', [ -916   484   422   423]);
hold on
imagesc(mydll)
axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
axis off
axis ij

hold on
plot(datarun.cones.centers(:,1)*9, datarun.cones.centers(:,2)*9, 'or', 'markersize', 4)

%% Voronoi +  cone marks


datarun = make_voronoi_masks(datarun);
masks = datarun.cones.mosaic.voronoi_masks;
excludes = [];
indexes = 1:length(masks);
indexes = setdiff(indexes, excludes);
cone_map = index_masks(masks, num2cell(indexes));

min_neighbor_dist = 1.5;
max_self_dist     = 4;
spaced = space_voronoi_masks(datarun, min_neighbor_dist, max_self_dist);
cone_map = index_masks(spaced, num2cell(indexes));

tmp = imresize(cone_map, 4.5, 'method','nearest' );

% final
ttt = tmp;
ttt(ttt>0) = 1;
figure
set(gcf, 'position', [ -916   484   422   423]);
hold on
colormap gray
imagesc(ttt)

axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
axis ij
axis off
plot(datarun.cones.centers(:,1)*9, datarun.cones.centers(:,2)*9, 'or', 'markersize', 4)



figure
imagesc(tmp)
axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])

% as specified by EJ

tt = zeros(2700, 2700, 3);
tt(1:2694,1:2694,:) =mydll/max(mydll(:))/1.1;

% tt(4:2697,4:2697,:) =mydll/max(mydll(:))/1.1;
ttt = tmp;
ttt(ttt>0) = 0.3;
tt(:,:,2) = tt(:,:,2)+ttt;
% figure
imagesc(tt)

axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
hold on
plot(datarun.cones.centers(:,1)*9, datarun.cones.centers(:,2)*9, 'or', 'markersize', 4)

%% 6: marked cones
cellIDs
ttt = tmp;
ttt(ttt>0) = 1;
figure
colormap gray
imagesc(ttt)
hold on
cnt = 1;
col = 'rgbkrgbkrgbkrg'
marks = 'ox+ox+ox+ox+ox+'
axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
for mc = cellIDs
    cw = datarun.cones.weights(:,mc);
    cw = cw/max(cw);
%     cw(cw<=0) =0.1;
    
    % figure
    % hold on
    % axis ij
    % for i=1:1032
    %     plot(datarun.cones.centers(i,1), datarun.cones.centers(i,2), 'or', 'markersize', cw(i)*5)
    % end
    
totcones = 0;
    for i=1:1032
        if cw(i)>0.25
            totcones = totcones+1;
            plot(datarun.cones.centers(i,1)*9, datarun.cones.centers(i,2)*9, [marks(cnt), col(cnt)], 'markersize', cw(i)*10)
        end
    end
    totcones
    cnt = cnt+1;
    drawnow
end


%% 6: colored cones
cellIDs
ttt = tmp;
ttt(ttt>0) = 1;
figure
set(gcf, 'position', [ -916   484   422   423]);
hold on
colormap gray
imagesc(ttt)
cnt = 1;
colors = [0.5 0.5 0.5; 1 0.5 0; 0.5 0 1;  1 1 0;0 1 0;  0 0 1;  1 0 0;1 0 1; 0 1 1; 0.6 0.65 0; ...
    0.8 0.8 0; 0.8 0 0.8; 0 0.8 0.8; 0 0.3 0.8; 1 1 0.3; 0.9 0.3 0.8; 0.1 0.1 0.3; 0.8 0.8 0.8; 0.1 0.1 0.1; 0 0.6 0; ];
axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
lay1 = tmp;
lay2 = tmp;
lay3 = tmp;
pk = cell(1,length(cellIDs));
for mc = cellIDs
    cw = datarun.cones.weights(:,mc);
    cw = cw/max(cw);
    newm = zeros(2700);
    totcones = 0;
    for i=1:1032
        if cw(i)>0.17
            ind = tmp(round(datarun.cones.centers(i,2)*9),round(datarun.cones.centers(i,1)*9));
            if ind>0
                totcones = totcones+1;                
                lay1(tmp==ind) = colors(totcones,1);
                lay2(tmp==ind) = colors(totcones,2);
                lay3(tmp==ind) = colors(totcones,3);
                plot(datarun.cones.centers(i,1)*9, datarun.cones.centers(i,2)*9, '.', 'color', colors(totcones,:), 'markersize', 30)
                newm(tmp==ind) = 1;
            end
        end
    end
    [r,c] = find(newm);
    a=convhull(r,c);
    pk{cnt} = [c(a),r(a)];
%     figure
    plot(c(a),r(a), 'linewidth',3)
%     totcones
    cnt = cnt+1;
    drawnow
end

axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
axis off



comb = zeros(2700,2700,3);
comb(:,:,1) = lay1;
comb(:,:,2) = lay2;
comb(:,:,3) = lay3;
figure
set(gcf, 'position', [ -916   484   422   423]);
hold on
imagesc(comb)
hold on
for i=1:length(cellIDs)
    plot(pk{i}(:,1), pk{i}(:,2), 'linewidth',3);
end
axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
axis ij
axis off

figure
set(gcf, 'position', [ -916   484   422   423]);
hold on
for i=1:length(cellIDs)
    plot(pk{i}(:,1), pk{i}(:,2), 'r','linewidth',3);
end
axis([130,210,70,155]*9)
set(gca,'dataaspectratio', [1 1 1])
axis ij
axis off
%%

% imagesc(mydll(90*9:170*9,138*9:201*9))
plot_rf_summaries(datarun, {cell_type}, 'scale', 9,'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'r')


figure
for i=1:14
    sta = squeeze(datarun.stas.stas{cellIDs(i)});
    subplot(4,4,i)
    colormap gray
    imagesc(sta(:,:,4))
    plot_rf_summaries(datarun, {cell_type},'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'b', 'contour_widths', 2)
    plot_rf_summaries(datarun, datarun.cell_types{cell_type}.cell_ids(i),'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'r')
    axis([135,205,80,160])
end

tmp = imresize(mydll,1/9, 'method', 'nearest');
figure
imagesc(tmp)

figure
for i=1:14
    sta = squeeze(datarun.stas.stas{cellIDs(i)});
    sta = -double(sta(:,:,4));
    sta = sta/max(abs(sta(:)));
%     sta = sta+0.5;
%     sta(sta<0.7) = 0;
    comb = zeros(300,300,3);
    comb(:,:,2) = tmp(:,:,1)/max(tmp(:));
    comb(:,:,3) = tmp(:,:,1)/max(tmp(:));
    comb(:,:,1) = sta;
    subplot(4,4,i)
%     colormap gray
    imagesc(comb)
    plot_rf_summaries(datarun, {13},'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'b')
    plot_rf_summaries(datarun, datarun.cell_types{13}.cell_ids(i),'clear', false, 'label', false, 'plot_fits', true, 'fit_color', 'r')
    axis([135,205,80,160])
end

