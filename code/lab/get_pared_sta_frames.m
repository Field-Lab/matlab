function datarun = get_pared_sta_frames(datarun, cell_spec, frames_to_keep)
% get_pared_sta_frames     Keep only a limited number of frames for each STA
%
% usage:  datarun = get_pared_sta_frames(datarun, cell_spec, frames_to_keep)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells
%      frames_to_keep - number of frames to keep
%
% outputs:    datarun - same struct as input, with modification to datarun.stas.stas{:}
%




% BODY OF THE FUNCTION

% get list of cells
cell_nums = get_cell_indices( datarun, cell_spec);

% go through through list of cells
for cell_num = cell_nums

    % get full length STA
    full_sta = datarun.stas.stas{cell_num};
    
    % skip if there is no STA
    if isempty(full_sta)
        continue
    end
    
    % keep only last frames
    datarun.stas.stas{cell_num} = datarun.stas.stas{cell_num}(:,:,:,end-frames_to_keep + 1:end);
    
end

