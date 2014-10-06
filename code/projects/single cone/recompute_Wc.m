function Wc = recompute_Wc(datarun,bcf,bcf_params,final_magic_number,extra_cones)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BAYESIAN CONE FINDING, ADDENDUM      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% recompute_Wc      write out new Wc
%
%
% usage:  Wc = recompute_Wc(datarun,bcf,bcf_params,final_magic_number,extra_cones)
%
% arguments:     datarun - datarun struct
%                    bcf - struct of results from bayesian cone finding
%             bcf_params - struct of parameters from bayesian cone finding
%     final_magic_number - the magic number to use
%            extra_cones - struct array, same format as kernel_spec in make_cone_weights_matrix
%
%
%
%
% 2010-08  gauthier
%




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
if 0

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




% SAVE OUT SOME OTHER STUFF

% cone weights matrix
%save([save_file_path 'Wc'],'Wc')


% note duration
fprintf('\n%recomputed Wc in %0.1f min\n\n\n',etime(clock,start_time_legitimize)/60)

