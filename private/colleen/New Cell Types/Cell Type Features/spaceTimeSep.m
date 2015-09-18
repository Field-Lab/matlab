% space time separability

clear all
% close all
clc
tic

dataparam.date='2015-08-17-1';
dataparam.concatname='d01-29-norefit';

% dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.file_name = [dataparam.date, '/',dataparam.concatname, '/data018/data018'];

% dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-20-1-0.48-22222.xml';

dataparam.cell_type = { 'ON midget', 'ON parasol', 'ON large 1', 'ON large 2'};
cmap = hsv(length(dataparam.cell_type));
fitparam.num_frames = 30;
        fitparam.pad_factor = 5;

%% END OF INPUT
dataparam.folder = dataparam.concatname;
% file path to save data and pictures
dataparam.filepath=['/Users/colleen/Desktop/SpaceTimeSep/',dataparam.date,'/',dataparam.concatname,'/data007/'];
if ~exist([dataparam.filepath],'dir')
    mkdir([dataparam.filepath]);
end

datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name, '.sta'];


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

for type = 1:size(dataparam.cell_type,2)
%     dataparam.folder = dataparam.cell_type{type};
    % file path to save data and pictures
%     dataparam.filepath=['/Users/colleen/Desktop/SpaceTimeSep/',dataparam.date,'/',dataparam.concatname,'/'];
%     if ~exist([dataparam.filepath],'dir')
%         mkdir([dataparam.filepath]);
%     end
    
    cell_type_index= zeros(1,size(dataparam.cell_type,2));
    for num_cell_types = 1:size(dataparam.cell_type,2)
        for i = 1:size(datarun.cell_types,2)
            right_cell_type = strcmpi(datarun.cell_types{i}.name, dataparam.cell_type{num_cell_types}); % case insensitive
            if right_cell_type == 1;
                cell_type_index(num_cell_types) = i;
                break
            end
            cell_type_index(num_cell_types) = 0;% couldn't find the right cell type
        end
        
    end
    
    % Set the cell_specification to all the cell of the inputted type
    dataparam.cell_specification = datarun.cell_types{cell_type_index(type)}.cell_ids;
    cell_indices = get_cell_indices(datarun, dataparam.cell_specification);
    cell_ids=datarun.cell_ids(cell_indices);
    
    for rgc = 1:length(cell_indices)
        
        sta = double(datarun.stas.stas{cell_indices(rgc)});
        if sum(isnan(sta(:))) ~= 0 || isnan(datarun.stas.fits{cell_indices(rgc)}.sd(1))
             parameter1(type,rgc) = nan;
            parameter2(type,rgc) = nan;
            parameter3(type,rgc) = nan;
                parameter4(type,rgc) = nan;
        else
        
        % figure; imagesc(sta(:,:,2,52));
        if size(sta,3)>1
            sta_one = squeeze(sta(:,:,2,:));
        else 
            sta_one = squeeze(sta);
        end
        
        
        
        scale = [1 1];
        aspect_ratio = 1;
        fit = datarun.stas.fits{cell_indices(rgc)};
        
        xdiff = fit.sd(1) / 2 * fitparam.pad_factor;
        ydiff = fit.sd(2) / 2 * fitparam.pad_factor;
        
        ydiff = max(ydiff, xdiff/aspect_ratio);
        xdiff = max(xdiff, ydiff*aspect_ratio);
        
        xstart = fit.mean(1) - xdiff;
        xend   = fit.mean(1) + xdiff;
        ystart = fit.mean(2) - ydiff;
        yend   = fit.mean(2) + ydiff;
        bounds = round([scale(1) scale(1) scale(2) scale(2)].*[xstart xend ystart yend]);

        if bounds(1)< 1
            bounds(1) = 1;
        end
         if bounds(3)< 1
            bounds(3) = 1;
         end
         if bounds(2)> size(sta_one,2)
            bounds(2) = size(sta_one,2);
         end
         if bounds(4)> size(sta_one,1)
            bounds(4) = size(sta_one,1);
        end
%             figure; 
%         imagesc(sta_one(bounds(3):bounds(4), bounds(1):bounds(2), 27))
        sta_cropped = sta_one(bounds(3):bounds(4), bounds(1):bounds(2), size(sta_one,3)-fitparam.num_frames+1:size(sta_one,3));
        sta_reshape = reshape(sta_cropped, size(sta_cropped,1)*size(sta_cropped,2),fitparam.num_frames);
    
        [U,D,V] = svd(sta_reshape);
        %     if type == 1 || type == 2
        %         figure; plot(diag(D), 'o')
        %         title(num2str(type))
        %         set(gca,'ylim', [0 0.75])
        %     end
        %     title([ 'cell id ' num2str(cell_ids(rgc))])
        if type ==1
        figure; plot(diag(D), 'o-')
        end
        parameter1(type,rgc) = D(1,1)/D(2,2);
        parameter2(type,rgc) = D(1,1)/D(3,3);
        parameter3(type,rgc) = D(1,1)/sum(D(:));
                parameter4(type,rgc) = (D(2,2)+D(1,1))/sum(D(:));
        end
        
% 
%                     eigenvalues = diag(D);
%         total_var = sum(eigenvalues);
%         for d = 1:length(eigenvalues)
%             var_explained(d) = sum(eigenvalues(1:d))/total_var;
%         end
%         parameter5{type,rgc} = var_explained;
        
    end
end

save([dataparam.filepath, '/', dataparam.folder],'parameter1', 'parameter2', 'parameter3', 'parameter4' )

toc


FigHandle = figure('Position', [100, 100, 1600, 1000]);
subplot(2,2,1)
for i = 1:size(dataparam.cell_type,2)
    temp = nonzeros(parameter1(i,:));
    plot(i*ones(size(temp,1)), temp,'.', 'MarkerSize', 20, 'Color', cmap(i,:));
    hold on
    set(gca, 'xlim', [0 size(dataparam.cell_type,2)+1]);
end
set(gca, 'xtick', 1:size(dataparam.cell_type,2));
set(gca, 'xticklabel', dataparam.cell_type, 'fontsize', 12)
title('\lambda_1 / \lambda_2')


subplot(2,2,2) 
for i = 1:size(dataparam.cell_type,2)
    temp = nonzeros(parameter2(i,:));
    plot(i*ones(size(temp,1)), temp, '.', 'MarkerSize', 20,'Color', cmap(i,:));
    hold on
    set(gca, 'xlim', [0 size(dataparam.cell_type,2)+1]);
end
set(gca, 'xtick', 1:size(dataparam.cell_type,2));
set(gca, 'xticklabel', dataparam.cell_type, 'fontsize', 12)
title('\lambda_1 / \lambda_3')



subplot(2,2,3)
for i = 1:size(dataparam.cell_type,2)
    temp = nonzeros(parameter3(i,:));
    plot(i*ones(size(temp,1)), temp, '.', 'MarkerSize', 20,'Color', cmap(i,:));
    hold on
    set(gca, 'xlim', [0 size(dataparam.cell_type,2)+1]);
end
set(gca, 'xtick', 1:size(dataparam.cell_type,2));
set(gca, 'xticklabel', dataparam.cell_type, 'fontsize', 12)
title('\lambda_1 / sum(\lambda)')


subplot(2,2,4)
for i = 1:size(dataparam.cell_type,2)
    temp = nonzeros(parameter4(i,:));
    plot(i*ones(size(temp,1)), temp, '.', 'MarkerSize', 20,'Color', cmap(i,:));
    hold on
    set(gca, 'xlim', [0 size(dataparam.cell_type,2)+1]);
end
set(gca, 'xtick', 1:size(dataparam.cell_type,2));
set(gca, 'xticklabel', dataparam.cell_type, 'fontsize', 12)
title('(\lambda_1 + \lambda_2)/ sum(\lambda)')


%     figure; 

%     
% cmap = hsv(size(dataparam.cell_type,2));
% for i = 1:size(dataparam.cell_type,2)
%     temp = find(~cellfun(@isempty,parameter5(i,:)));
%     for j = 1:length(temp)
%         semilogy(parameter5{i,j}, 'color',cmap(i,:))
%         hold on
%     end
%     
% end


suptitle({dataparam.date; dataparam.concatname})
    export_fig([dataparam.filepath, 'SpaceTimeSep'], '-pdf')

% print(gcf,'-dpdf',sprintf('%s%s%s.pdf',[dataparam.filepath,'/'],['SpaceTimeSep']));

