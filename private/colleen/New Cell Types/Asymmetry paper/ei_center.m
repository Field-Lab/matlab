function spacing = ei_center(datarun, cell_ids)
% Position of center of EI

% dataparam.date='2006-05-05-0';
% dataparam.concatname='data000';
% 
% 
% % dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
% dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
% dataparam.save_path = ['/Volumes/Lab/Users/crhoades/JitterMovie/', dataparam.date, '/', dataparam.concatname, '/'];
% % list specific cell (1), or run for a whole cell type (0)
% select_cells = 0;
% if select_cells == 1
%     dataparam.vision_ids= [2] %ON parasol
% end
% dataparam.cell_specification = {'ON parasol'};
% 
% 
% %% END OF INPUT
% % dataparam.folder = dataparam.cell_type{1};
% % file path to save data and pictures
% % dataparam.filepath=['/Users/colleen/Desktop/Fitting/',dataparam.date,'/',dataparam.concatname,'/data023/'];
% % if ~exist([dataparam.filepath,dataparam.folder],'dir')
% %     mkdir([dataparam.filepath,dataparam.folder]);
% % end
% 
% % Right Movie
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name, '.neurons'];
% datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name, '.params'];
% datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name, '.sta'];
% datarun.names.rrs_ei_path = ['/Volumes/Analysis/', dataparam.file_name, '.ei'];
% 
% %% Load Data
% 
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
% opt.load_sta_params.save_rf = 1;
% % opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
% datarun=load_data(datarun,opt);
% 
% dataparam.cell_ids = get_cell_indices(datarun, dataparam.cell_specification);
% dataparam.vision_ids = get_cell_ids(datarun, dataparam.cell_specification);
% 
% % datarun.piece.array_id = 1505;
% datarun = load_ei(datarun, dataparam.vision_ids);
% % if size(datarun.ei.eis{1},1) ==519
% %     datarun = load_ei(datarun, dataparam.vision_ids, 'array_id', 1505);
% % end
% 
% 
% % dataparam.cell_ids = get_cell_indices(datarun,dataparam.cell_specification);
if size(datarun.ei.eis{1},1) == 519
    
    [xc,yc] = getElectrodeCoords519();
else
    [xc,yc] = getElectrodeCoords512();

end

if size(xc,1) == 1
    xc = xc';
    yc= yc';
end

for i = 1:length(cell_ids);
ei = datarun.ei.eis{cell_ids(i)};

% convert disk size from units of median nearest neighbor spacing to absolute units
params.max_scale = infer_electrode_spacing([xc, yc])/2;

%normalize
ei = ei/abs(min(ei(:)));


ei_frame = get_ei_max_frame(ei, 1, size(ei,2), 0);


waveforms_ordered = sort(ei_frame);
smallest = waveforms_ordered(3);
indices = find(ei_frame<=smallest);



% circle samples
if size(datarun.ei.eis{1},1) == 519
    circle_samples = [0:(2*pi/127):2*pi 2*pi];

else
    circle_samples = [0:(2*pi/125):2*pi 2*pi];

end

x_circle = cos(circle_samples);
y_circle = sin(circle_samples);



totalmass = 0;
totalx = 0;
totaly = 0;
for j = 1:3
    totalmass = totalmass + ei_frame(indices(j));
    totalx = totalx+ xc(indices(j),1)*ei_frame(indices(j));
    totaly = totaly+ yc(indices(j),1)*ei_frame(indices(j));
end
center(i,:) = [totalx/totalmass totaly/totalmass];


% draw a disk for each electrode not equal to 0
% figure
% scatter(xc,yc, 1*ones(size(xc)), 'k','filled')
% 
% hold on
% 
% for e = find(ei_frame ~= 0)'
%     
%     % size of circle
%     sd = params.max_scale * abs(ei_frame(e));
%     if sd > params.max_scale,
%         sd = params.max_scale;
%     end
%     
% 
%     X = sd * [x_circle; y_circle];
%     
%     % draw patch, or plot outline
%     
%     elec_handle = patch(X(1,:) +xc(e), X(2,:) +yc(e), [0 0 0], 'Facealpha', 0.5);
%     
%     
%     zdata = ones(size(get(elec_handle, 'XData')));
%     set(elec_handle, 'ZData', zdata);
%     
%     
%     elec_handles(e) = elec_handle;
% end
% 
% scatter(center(i,1), center(i,2), 300, 'r', 'filled')
% hold on;
% box on;
% axis equal;
% range=[min(xc) max(xc) min(yc) max(yc)]*1.1;
% 
% axis(range);
end

% figure; scatter(center(:,1), center(:,2))

[IDX, D] = knnsearch(center, center, 'K', 5);
D = D(:,2:end);
total = zeros(size(D,1),1);

for i = 1:size(D,1)
    if D(i,1) > 10
        total(i) = D(i,1);
    elseif D(i,2) > 10
        total(i) = D(i,2);
    elseif D(i,3) > 10
        total(i) = D(i,3);
    else
        total(i) = D(i,4);
    end
end

[N, bin_centers] = hist(total,7);

[val, ind] = max(N);
spacing = bin_centers(ind);




