function datarun = copy_images_to_stacks(datarun, varargin)
% COPY_IMAGES_TO_STACKS     Copy image fields read from index file to new stacks system
% usage: datarun = copy_images_to_stacks(datarun, opts)
%
% opts:     basepath    server_data_path    Path to images; as of 2010-09 this is jacob
%
% This is run automatically in load_index
%
% 2010-09 phli
%

if ~isfield(datarun, 'piece')
    return;
end

opts = inputParser;
opts.addParamValue('basepath', server_data_path);
opts.parse(varargin{:});
opts = opts.Results;

piece_info = datarun.piece;
% [stack_cell_labels, stack_row_labels] = stack_labels();
% pm_row = stack_row_labels.photographic_mapping;

mapping_names = [{'array_edges' 'array'} load_alignment_image()];
field_names = collect(mapping_names, @(nam) ['photographic_mapping_' nam]);
field_names = [field_names {'alive_montage_lores' 'alive_montage_hires'}];

for i = 1:length(field_names)
    field_name = field_names{i};

    % Do not overwrite
    if ~stackempty(datarun, field_name)
        continue
    end
    
    if isfield(piece_info, field_name)
        stack = build_stack(piece_info.(field_name), 'basepath', opts.basepath, 'name', field_name);
        if ~isempty(stack)
            datarun = set_stack(datarun, stack);
        end
    end
end