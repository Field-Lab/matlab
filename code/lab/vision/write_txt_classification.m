function write_txt_classification(classification, savepath, varargin)
% WRITE_TXT_CLASSIFICATION      Write Vision format classification.txt file from Matlab classification struct
% usage: write_txt_classification(classification, savepath)
%
% CLASSIFICATION is a nested struct with fields: "name", "subclasses",
% "cells".  NAME is the string name of the cell class.  SUBCLASSES is a
% substruct with the same format, or an empty struct.  CELLS is an integer
% vector with the cell_ids of the cells in this class.
%
% To be compatible with vision, the top category is
% conventionally called "All".
%
% SAVEPATH is the full path and filename of the classification.txt file to
% output.
%
% OPTS:
%   verbose     false
%
% See also RESHAPE_TXT_CLASSIFICATION
%
% 2011-05 phli
%

opts = inputParser;
opts.addParamValue('verbose', false);
opts.parse(varargin{:});
opts = opts.Results;


% Show overwrite existing txt dialog
if exist(savepath, 'file')
    overwrite = questdlg(['File ' savepath ' exists; overwrite?'], 'Overwrite existing classification?', 'Okay', 'Cancel', 'Cancel');
    if strcmp(overwrite, 'Cancel')
        return
    end
end


reshaped = reshape_txt_classification(classification);


fid = fopen(savepath, 'w');
if fid < 0
    error(['Failed to open file ' savepath ' for writing!']);
end


try
    disp(['Saving classification to: ' savepath '...']);
    for i = 1:length(reshaped)
        if isempty(reshaped{i})
            continue;
        end
        
        out = sprintf('%d  %s', i, reshaped{i});
        if opts.verbose
            disp(out);
        end
        fprintf(fid, [out '\n']);
    end
catch e
    % Ensure
    fclose(fid);
    rethrow(e);
end


fclose(fid);