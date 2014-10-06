function ei = get_ei(datarun, cell_id)
% get_ei     return an EI
%
% usage:  ei = get_ei(datarun, cell_id)
%
% arguments:     datarun - datarun struct
%                cell_id - cell id
%
% outputs:     ei - ExF matrix.  E = # electrodes, F = # frames
%
%
% 2010-02 gauthier
%


% get cell index
cell_index = get_cell_indices(datarun,cell_id);

% get the EI
ei = datarun.ei.eis{cell_index};


% if it's not loaded, load it
if ~isempty(ei)
    return
else
    % load EI data
    if isfield(datarun,'ei') && isfield(datarun.ei,'java_ei') && ~isempty(datarun.ei.java_ei)
        ei = datarun.ei.java_ei.getImage(datarun.cell_ids(cell_index));
    else
        % NOTE: the java ei object should only be loaded once. this is because matlab "opens a file" each time
        % a java object is loaded, and it can only open a limited number of files.  note that these "open files" for java objects
        % are DIFFERENT than matlab's standard list of open files.  as far as JG can tell, you can not check how many
        % open java files there are, nor is the maximum number documented.
        error(sprintf(['java object of the ei file must be loaded first.  Use this command:\n'...
            '   datarun = load_ei(datarun,[]);']))
    end
    
    % throw out the error
    ei = squeeze(ei(1,:,:));
    
    % throw out the triggers
    ei = ei(2:end,:);
end