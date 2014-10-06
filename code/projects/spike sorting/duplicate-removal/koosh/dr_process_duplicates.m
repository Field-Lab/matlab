% this algorithm compares spike trains and merges them if they appear to be
% duplicates of one another. see the paper in the /references folder.
%
% this function isn't done, but it's a good starting point for further
% exploration.
%
%
function processed_cells = dr_process_duplicates(datarun, samplingRate)

contam_thresh  = 0.0;
compare_thresh = 0.25;

% find all cells with contamination <= threshold
% (assumed as a precondition before calling this function)

% populate unprocessed_cells queue with spike trains
unprocessed_cells.s = cell(length(datarun.spikes),1);
for ii = 1:length(datarun.spikes)
    unprocessed_cells.s{ii} = datarun.spikes{ii} * samplingRate;
    unprocessed_cells.e(ii) = datarun.channels(ii);
end

% initialize processed_cells to be empty
processed_cells.s = {};
processed_cells.e = [];

% initialize killed_cells to be empty
killed_cells.s = {};
killed_cells.e = [];

% run koosh
while ~isempty(unprocessed_cells.s)
    
    % get first item from the queue
    disp(sprintf('unprocessed length: %d', length(unprocessed_cells.s)))
    current.s  = unprocessed_cells.s{1};
    current.e  = unprocessed_cells.e(1);
    unprocessed_cells.s(1) = [];
    unprocessed_cells.e(1) = [];
    
    % flag if we've merged at least once this round
    mergeFlag = false;
    
    % compare current cell against processed cells
    for pair_index = 1:length(processed_cells.s)
        
        % we are removing items from processed cells
        if pair_index > length(processed_cells.s)
            break;
        end
        
        paired.s = processed_cells.s{pair_index};
        paired.e = processed_cells.e(pair_index);
        
        % merge combines cells if they are sufficiently similar
        s.a = current.s;
        s.b = paired.s;
        e.a = current.e;
        e.b = paired.e;
        if isempty(s.a) || isempty(s.b), continue; end
        [merged] = dr_compare_spike_trains(s, e, datarun.duration, compare_thresh, samplingRate);
        
        % merged is empty if the comparison function has indicated that the
        % two spike trains are not similar at the threshold specified
        if isempty(merged), continue; end
        
        disp(sprintf('%d and %d have been merged!', length(s.a), length(s.b)))
        disp(sprintf('%d spikes gained!', length(merged.s)-max(length(s.a),length(s.b))))
        mergeFlag = true;
        processed = false;
        
        % compare merged cell to current cell
        if isequal(merged.s, current.s)
            processed = true;
        end

        % compare merged cell to unprocessed cell list
        if processed == false
            for u_spikes = 1:length(unprocessed_cells.s)
                u = unprocessed_cells.s{u_spikes};
                if length(merged.s) ~= length(u)
                    continue;
                end
                if isequal(merged.s, u)
                    processed = true;
                    break;
                end
            end
        end

        % compare merged cell to processesed list
        if processed == false
            for p_spikes = 1:length(processed_cells.s)
                p = processed_cells.s{p_spikes};
                if length(merged.s) ~= length(p)
                    continue;
                end
                if isequal(merged.s, p)
                    processed = true;
                    break;
                end
            end
        end

        % if merged cell is in none of those lists, add it to unmerged list
        % for reprocessing
        if processed == false
            % add merged cell to the queue
            unprocessed_cells.s{end + 1} = merged.s; %#ok<AGROW>
            unprocessed_cells.e(end + 1) = merged.e; %#ok<AGROW>
            
            % remove children from the processed list, put it in killed
            % list so it can be checked again. just because one possible
            % merge has occured doesn't mean other valid merges can't occur
            % since we are not necessarily dealing with equivalence classes
            killed_cells.s{end+1} = paired.s; %#ok<AGROW>
            killed_cells.e(end+1) = paired.e; %#ok<AGROW>
            processed_cells.s(pair_index) = []; %#ok<AGROW>
            processed_cells.e(pair_index) = []; %#ok<AGROW>
            pair_index = pair_index - 1; %#ok<FXSET>
        end
    end
    
    % compare current cell against killed cells
    for pair_index = 1:length(killed_cells.s)
        
        % we are removing items from killed cells
        if pair_index > length(killed_cells.s)
            break;
        end
        
        paired.s = killed_cells.s{pair_index};
        paired.e = killed_cells.e(pair_index);
        
        % merge combines cells if they are sufficiently similar
        s.a = current.s;
        s.b = paired.s;
        e.a = current.e;
        e.b = paired.e;
        if isempty(s.a) || isempty(s.b), continue; end
        [merged2] = dr_compare_spike_trains(s, e, datarun.duration, compare_thresh, samplingRate);
        
        % merged is empty if the comparison function has indicated that the
        % two spike trains are not similar at the threshold specified
        if isempty(merged2) || isequal(merged, merged2)
            continue
        end
        
        disp(sprintf('%d and %d have been merged!', length(s.a), length(s.b)))
        disp(sprintf('%d spikes gained!', length(merged2.s)-max(length(s.a),length(s.b))))
        processed = false;
        
        % compare merged cell to current cell
        if isequal(merged2.s, current.s)
            processed = true;
        end

        % compare merged cell to unprocessed cell list
        if processed == false
            for u_spikes = 1:length(unprocessed_cells.s)
                u = unprocessed_cells.s{u_spikes};
                if length(merged2.s) ~= length(u)
                    continue;
                end
                if isequal(merged2.s, u)
                    processed = true;
                    break;
                end
            end
        end

        % compare merged cell to processesed list
        if processed == false
            for p_spikes = 1:length(processed_cells.s)
                p = processed_cells.s{p_spikes};
                if length(merged2.s) ~= length(p)
                    continue;
                end
                if isequal(merged2.s, p)
                    processed = true;
                    break;
                end
            end
        end

        % if merged cell is in none of those lists, add it to unmerged list
        % for reprocessing
        if processed == false
            % add merged cell to the queue
            unprocessed_cells.s{end + 1} = merged2.s; %#ok<AGROW>
            unprocessed_cells.e(end + 1) = merged2.e; %#ok<AGROW>
        end
    end
    
    
    % add current cell to clean cell list if we didn't merge this round
    % otherwise add it to the killed list
    disp(sprintf('processed length: %d', length(processed_cells.s)))
    if mergeFlag == false
        processed_cells.s{end+1} = current.s; %#ok<AGROW>
        processed_cells.e(end+1) = current.e; %#ok<AGROW>
    else
        killed_cells.s{end+1} = current.s; %#ok<AGROW>
        killed_cells.e(end+1) = current.e; %#ok<AGROW>
    end
end

% phase 2 is greedy duplicate group selection

% or if merge, remove children (since we're getting larger)
% but this makes order matter

end