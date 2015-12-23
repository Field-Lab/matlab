function reshaped = reshape_txt_classification(classification, reshaped, parent_string)
% RESHAPE_TXT_CLASSIFICATION    Restructure Matlab cell classification nested structs into format ready to write to Vision classification.txt
% usage: reshaped = reshape_txt_classification(classification)
%
% CLASSIFICATION is a nested struct with fields: "name", "subclasses",
% "cells".  NAME is the string name of the cell class.  SUBCLASSES is a
% substruct with the same format, or an empty struct.  CELLS is an integer
% vector with the cell_ids of the cells in this class.
%
% This reshapes CLASSIFICATION into a format ready to write to a
% classification.txt file.  Output RESHAPED is a cell array with empty
% elements except for indices that are cell ids, where the elements will be
% the classification string, e.g. "All/ON/Midget".
%
% Calls itself recursively with extra arguments descend the nested
% classification structs.
%
% 2011-05 phli
%


% Extra arguments used for recursive calls
if nargin < 3
    parent_string = '';
end
if nargin < 2
    reshaped = {};
end


classification_string = [parent_string classification.name];

% Add cells
for i = 1:length(classification.cells)
    cell_id = classification.cells(i);
    reshaped{cell_id} = classification_string;
end

% Recursively add subclasses
for i = 1:length(classification.subclasses)
    reshaped = reshape_txt_classification(classification.subclasses(i), reshaped, [classification_string '/']);
end