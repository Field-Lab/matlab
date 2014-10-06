function datarun = merge_stacks(datarun, loaded)
% MERGE_STACKS      Combine loaded stacks data into datarun
% usage: datarun = merge_stacks(datarun, loaded)
%
% This is a comewhat complicated process.  The stacks array has a bunch of
% different stack objects, which are images or sets of images with common
% coordinates.  When stacks is loaded from disk and merged with datarun,
% there is the danger that the two may not be completely consistent.  So we
% first check through all the stacks to see where there might be
% conflicting entries in the 2D arrays.
%
% If inconsistencies are found, we try to handle them as nicely as
% possible.  Stacks that are inconsistent are not merged into the datarun.
% Furthermore, there are transform fields for each stack, that hold
% transforms in relation to every other stack.  So if one of the loaded
% stacks is inconsistent and not merged in, then we should not merge in the
% transforms associated with that stack for any of the other loaded stacks.
%
% 2010-10 phli
%


% First decide if the loaded stacks are consistent with datarun.stacks, and which ones might be inconsistent
[matched_all, matched_indices] = check_stacks_consistency(datarun.stacks, loaded.stacks);
if ~matched_all
    warning('Loaded stacks were not all consistent with datarun.stacks.  Loading selectively.');
end

% Now go through all the stacks and merge them together, skipping those
% that were found to be inconsistent
datarun_stacks_size = size(datarun.stacks);
loaded_stacks_size = size(loaded.stacks);
for i = 1:size(matched_indices, 1)
    for j = 1:size(matched_indices, 2)
        if ~matched_indices(i,j)
            warning(['Skipping inconsistent stack: ' num2str(i) ',' num2str(j)]);
            continue
        end
        
        % If loaded is blank we are done
        if any([i j] > loaded_stacks_size) || isempty(loaded.stacks{i,j})
            continue
        end
        
        % If datarun is blank, automatically use loaded
        if any([i j] > datarun_stacks_size) || isempty(datarun.stacks{i,j})
            datarun.stacks{i,j} = loaded.stacks{i,j};
            continue
        end

        % If neither is blank, merge them together, being careful to throw
        % out transforms that correspond to stacks that were found to be
        % inconsistent
        dstack = datarun.stacks{i,j};
        lstack = loaded.stacks{i,j};
        lstack = strip_inconsistent_transforms(lstack, matched_indices);
        datarun.stacks{i,j} = merge_stack(dstack, lstack);
    end
end



function [bool, indices] = check_stacks_consistency(stacks1, stacks2)
size1 = size(stacks1);
size2 = size(stacks2);
minsize = min([size1;size2]);
maxsize = max([size1;size2]);

% If stack only exists in one or the other, it is automatically consistent
indices = true(maxsize);

% Go through the stacks that exist in both and check for consistency
for i = 1:minsize(1)
    for j = 1:minsize(2)
        indices(i,j) = check_stack_consistency(stacks1{i,j}, stacks2{i,j});
    end
end

% If any index is false, return false overall
bool = all(indices(:));



function bool = check_stack_consistency(stack1, stack2)
% If one is empty, then they are consistent
if isempty(stack1) || isempty(stack2)
    bool = true;
    return
end

% If names exist but don't match, then not consistent
if ~isempty(stack1.name) && ~isempty(stack2.name) && ~strcmp(stack1.name, stack2.name)
    bool = false;
    return
end

% If paths match, then consistent
if stack1.paths == stack2.paths
    bool = true;
    return
end

% If data match and isn't blank then consistent
if ~all(stack1.data, @isempty) && stack1.data == stack2.data
    bool = true;
    return
end

% Otherwise, inconsistent
bool = false;



function stack = strip_inconsistent_transforms(stack, indices)
% The transform fields are cell arrays that have elements corresponding to
% other stacks.  If another stack was found to be inconsistent, then we
% don't want to keep that transform, so clear our selected elements of each
% transform field.
tform_fields = {'input_points' 'base_points' 'tforms' 'tforms_inv' 'tform_routes' 'tforms_approx'};
for i = 1:size(indices, 1)
    for j = 1:size(indices, 2)
        % If is was consistent, move on
        if indices(i,j)
            continue
        end
        
        % If not consistent, blank corresponding element of every transform field
        for k = 1:length(tform_fields)
            tform_field = tform_fields{k};
            if isfield(stack, tform_field)
                stack.(tform_field){i,j} = [];
            end
        end
    end
end



function dstack = merge_stack(dstack, lstack)
% Go through each field name and merge in loaded values if dstack was missing them.  But do not overwrite!
field_names = union(fieldnames(dstack), fieldnames(lstack));
for i = 1:length(field_names)
    field_name = field_names{i};

    % If dstack doesn't exist or is empty, go ahead and fill in with lstack
    if ~isfield(dstack, field_name) || (isempty(dstack.(field_name)) && ~isempty(lstack.(field_name)))
        dstack.(field_name) = lstack.(field_name);
    end

    % For the fields that are cell arrays, actually go through the whole
    % cell array and fill in missing values in dstack from lstack
    if iscell(lstack.(field_name)) && ~all(lstack.(field_name), @isempty)
        for j = 1:size(lstack.(field_name), 1)
            for k = 1:size(lstack.(field_name), 2)
                if ~isempty(lstack.(field_name){j,k}) && isempty(dstack.(field_name){j,k})
                    dstack.(field_name){j,k} = lstack.(field_name){j,k};
                end
            end
        end
    end
end