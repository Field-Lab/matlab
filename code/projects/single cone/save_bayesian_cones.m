function save_bayesian_cones(datarun,bcf,bcf_params,final_magic_number,save_suffix,old_cone_finding_flag,extra_cones, varargin)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BAYESIAN CONE FINDING STEP 5 OF 5    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% save_bayesian_cones     write out cones to disk
%
%
% usage:  save_bayesian_cones(datarun,bcf, final_magic_number)
%
% arguments:     datarun - datarun struct
%                    bcf - struct of results from bayesian cone finding
%     final_magic_number - the magic number to use
%            save_suffix - text appended to folder name where files are saved
%  old_cone_finding_flag - also plot cones from old cone finding?
%            extra_cones - struct array, same format as kernel_spec in make_cone_weights_matrix
%
%
%
%
% 2010-02  gauthier
%


opts = inputParser();
opts.addParamValue('fit_foa', 4);
opts.addParamValue('robust_std_method', 3);
opts.parse(varargin{:});
opts = opts.Results;



% use bayesian cone list to compute fits to each RF.  export text files

start_time_legitimize = clock;


% EXPAND PARAMETERS

fn = fieldnames(bcf_params);
for ff=1:length(fn)
    eval(sprintf('%s = bcf_params.%s;',fn{ff},fn{ff}))
end

fn = fieldnames(bcf);
for ff=1:length(fn)
    eval(sprintf('%s = bcf.%s;',fn{ff},fn{ff}))
end



% IDENTIFY CONES TO KEEP

switch 1
    case 1 % all cone types have the same magic number
        % commented out: GDF to allow S cones to have a different magic number
        keep_indices = (all_added_cones(:,7) + final_magic_number*all_added_cones(:,8)) > 0;

    case 2 % different magic number for different cones types
        
        %% INTRODUCED BY GDF TO ALLOW S CONES TO HAVE A DIFFERENT THRESHOLD
        % get indices to LM versus S cones
        LM_indices = (all_added_cones(:,3)) ~= 3;
        S_indices = (all_added_cones(:,3)) == 3;

        %         % get kept LM cones
        %         keep_indices = (all_added_cones(:,7) + final_magic_number(1)*all_added_cones(:,8)) > 0;
        %         LM_keep_indices = keep_indices & LM_indices;
        %
        %         % get kept S cones
        %         clear keep_indices
        %         keep_indices = (all_added_cones(:,7) + (final_magic_number(2))*all_added_cones(:,8)) > 0;
        %         S_keep_indices = keep_indices & S_indices;
        %
        %         % join indices for LM and S cones
        %         keep_indices = LM_keep_indices | S_keep_indices;
        
        fmns = zeros(size(all_added_cones,1),1);
        fmns(LM_indices) = final_magic_number(1);
        fmns(S_indices) = final_magic_number(2);
        
        keep_indices = (all_added_cones(:,7) + fmns.*all_added_cones(:,8)) > 0;
        
    otherwise
        error('')
        
end

% discard unkept cones
kept_cones = all_added_cones(keep_indices,:);

% note cone count
num_cones = size(kept_cones,1);

%  note centers
cone_centers = kept_cones(:,[4 5]);

% type
kernel_color_names = fieldnames(kernel_colors);
cone_types = char(num_cones,1);
for cc = 1:size(kernel_plot_colors,1);
    cone_types(kept_cones(:,3)==cc) = kernel_color_names{cc};
end 

% likelihood
cone_delta_likelihoods = kept_cones(:,7);

% show how many were kept
if isscalar(final_magic_number)
    fprintf('Kept %d of %d cones (magic number = %0.2f)\n',sum(keep_indices),size(all_added_cones,1),final_magic_number)
else
    fprintf('Kept %d of %d cones (magic numbers, LM = %0.2f, S = %0.2f)\n',...
        sum(keep_indices),size(all_added_cones,1),final_magic_number(1),final_magic_number(2))
end





% MAKE PLOT OF CONES

% likelihood surface
figure(dll_fig);clf;

% dll plus cones
imagesc(matrix_scaled_up((norm_image(   dll   )-0.5).^0.4,3),...
    'xdata',[1 datarun.stimulus.field_width]-(1-kernel_spacing),'ydata',[1 datarun.stimulus.field_height]-(1-kernel_spacing))
axis image;hold on

% old cones
% load up old cones
if old_cone_finding_flag
    if ischar(old_cone_finding_flag)
        datarun = import_single_cone_data(datarun,old_cone_finding_flag);
    else
        datarun = import_single_cone_data(datarun,datarun.names.nickname);
    end
    plot_cone_mosaic(datarun,'fig_or_axes',gca,'clear',0,'bg_color',[],'drawnow',true,'cone_size',6,'type_colors',...
        [eye(3); 0 0 0] * .3 + 0.2)
end

% new cones
for cc=1:size(kernel_plot_colors,1)
    cns = (kept_cones(:,3) == cc);
    plot(kept_cones(cns,4),kept_cones(cns,5),'.','Color',kernel_plot_colors(cc,:),'MarkerSize',3)
end




% MAKE CONE WEIGHTS MATRIX
 
% size
rf_size = [datarun.stimulus.field_height datarun.stimulus.field_width];

% get locations
kernel_spec = struct('center',mat2cell(cone_centers,repmat(1,sum(keep_indices),1)));

% cone type
for cc = 1:size(kernel_plot_colors,1); [kernel_spec(kept_cones(:,3)==cc).type] = deal(kernel_color_names{cc}); end %#ok<USENS>

% cone radius
for cc = 1:size(kernel_plot_colors,1); [kernel_spec(kept_cones(:,3)==cc).radius] = deal(kernel_radii(cc)); end

% add in extra cones
if exist('extra_cones','var')
    kernel_spec = [kernel_spec extra_cones];
end

% make matrix
Wc = make_cone_weights_matrix(rf_size,kernel_spec,kernel_colors);




% GET CONE WEIGHTS IN EACH RGC


% get list of cell indices
cell_indices = get_cell_indices(datarun,regression_cell_spec);

% note when it started
fprintf('\nComputing cone weights in %d RFs',length(cell_indices));
start_time_regress = clock; 

% initialize
cone_weights = zeros(num_cones,length(cell_indices));

% go through list of cells
for cc = 1:length(cell_indices)
    
    fprintf('.')

    % get summary frame
%     temp_frame = length(datarun.stas.time_courses{1}(:,1)) -1;
%     [rf, junk] = get_sta(datarun,datarun.cell_ids(cell_indices(cc)), 'frames', temp_frame);
     rf = get_rf(datarun,datarun.cell_ids(cell_indices(cc)));

    if isempty(rf)
        continue
    end

    % reshape for the regression
    rf = reshape(rf,[],1);

    % put in units of SNR
    rf = rf / robust_std(rf, opts.robust_std_method);
    
    if size(rf,1)<size(Wc,1) 
        rf=repmat(rf,3,1);
    end

    % regress to get the cone weights

    cone_weights(:,cc) = Wc\rf; 
    %for cn = 1:size(Wc,2)
    %    cone_weights(cn,cc) = full(Wc(:,cn))\rf;
    %end
end

% display how long it took
fprintf('\n   done (%0.1f seconds)\n',etime(clock,start_time_regress));
    
% note list of cell ids
regression_cell_ids = datarun.cell_ids(cell_indices);

% clean up 
clear start_time_regress rf cc cell_indices





% FIT EACH RGC

clear cone_info
cone_info.cone_weights = cone_weights;
cone_info.cone_centers = cone_centers;
cone_info.cell_ids = get_cell_ids(datarun,regression_cell_spec);

rf_cone_fits = fit_cone_rfs(datarun,regression_cell_spec,'cone_info',cone_info,'foa_profile',[],'fit_radius',150,'verbose',1, 'foa_2d', opts.fit_foa);



% SAVE OUT CONE INFO


% load data into a struct
all_data = struct;
all_data.cone_weights = cone_weights;
all_data.cone_types = cone_types;
all_data.cone_centers = cone_centers;
all_data.rf_cone_fits = rf_cone_fits;
all_data.cone_likelihoods = cone_delta_likelihoods;

% meaningless data that must be loaded to have the correct number of fields
all_data.cone_rgb = ones(num_cones,3);
all_data.cone_types_kmeans = cone_types;
all_data.cone_types_em = cone_types;

path2data=datarun.names.rrs_prefix;
tmp=regexp(path2data,'data');
path2data=path2data(1:tmp(1)-1);

% save most info to text files
if isscalar(final_magic_number)
    save_file_path = [path2data datarun.names.short_name sprintf('-bayes-msf_%0.2f-%s/',final_magic_number,save_suffix)];
else
    save_file_path = [path2data datarun.names.short_name sprintf('-bayes-msf_LM_%0.2f-msf_S_%0.2f-%s/',...
        final_magic_number(1),final_magic_number(2),save_suffix)];
end
mkdir(save_file_path)
export_single_cone_data(datarun,regression_cell_spec,all_data,save_file_path)
clear all_data

% load data into datarun (needed for plotting stuff below)
datarun = import_single_cone_data(datarun,save_file_path);





% SAVE OUT SOME OTHER STUFF

% cone weights matrix
save([save_file_path 'Wc'],'Wc')

% all cones
save([save_file_path 'results'],'bcf')

% analysis parameters
% all_params = struct;
% all_params.cone_finding_cell_spec = cone_finding_cell_spec;
% all_params.magic_number = magic_number;
% all_params.final_magic_number = final_magic_number;
% all_params.q = q;
% all_params.kernel_spacing = kernel_spacing;
% all_params.kernel_radii = kernel_radii;
% all_params.LM_MM = LM_MM;
% all_params.LM_S = LM_S;
% all_params.S_S = S_S;
% all_params.relevant_region_x = relevant_region_x;
% all_params.relevant_region_y = relevant_region_y;
% all_params.num_iter = num_iter;
% all_params.regression_cell_spec = regression_cell_spec;
% all_params.num_iter = rel_radius; % added 2010-01
save([save_file_path 'parameters'],'bcf_params')





% MAKE PLOTS TO SAVE

% likelihood surface
figure(dll_fig);clf;

% make correct units
set(gcf,'PaperUnits','centimeters')
xPaper = 30; yPaper = 30;
xSize = xPaper; ySize = yPaper;
xLeft = (xPaper-xSize)/2; yTop = (yPaper-ySize)/2;
set(gcf,'PaperSize',[xPaper yPaper])
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])

% dll only

mydll=matrix_scaled_up((norm_image(   dll   )-0.5).^0.4,3);

save([save_file_path 'mydll'],'mydll')


subplot(211);
imagesc(matrix_scaled_up((norm_image(   dll   )-0.5).^0.4,3),...
    'xdata',[1 datarun.stimulus.field_width]-(1-kernel_spacing),'ydata',[1 datarun.stimulus.field_height]-(1-kernel_spacing))
axis image;hold on

% dll plus cones
subplot(212);
imagesc(matrix_scaled_up((norm_image(   dll   )-0.5).^0.4,3),...
    'xdata',[1 datarun.stimulus.field_width]-(1-kernel_spacing),'ydata',[1 datarun.stimulus.field_height]-(1-kernel_spacing))
axis image;hold on
 
% old cones
% load up old cones
if old_cone_finding_flag
    if ischar(old_cone_finding_flag)
        datarun = import_single_cone_data(datarun,old_cone_finding_flag);
    else
        datarun = import_single_cone_data(datarun,datarun.names.nickname);
    end
    plot_cone_mosaic(datarun,'fig_or_axes',gca,'clear',0,'bg_color',[],'drawnow',true,'cone_size',6,'type_colors',...
        [eye(3); 0 0 0] * .3 + 0.2)
end

% new cones
for cc=1:size(kernel_plot_colors,1)
    cns = (kept_cones(:,3) == cc);
    plot(kept_cones(cns,4),kept_cones(cns,5),'.','Color',kernel_plot_colors(cc,:),'MarkerSize',3)
end

% print
print(dll_fig,[save_file_path 'dll_surface'],'-dpdf')




% note duration
fprintf('\n%s: legitimized and saved everything in %0.1f min\n\n\n',datarun.names.nickname,etime(clock,start_time_legitimize)/60)

