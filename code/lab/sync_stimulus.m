function datarun = sync_stimulus(datarun, stimulus, source)
% SYNC_STIMULUS     Ensure that datarun has consistent stimulus information
%
% if datarun.stimulus.<field> doesn't match the provided stimulus, give an error
% if datarun.stimulus.<field> doesn't exist, set it to match the provided stimulus
%
%
% usage:  datarun = sync_stimulus(datarun, stimulus, source)
%
%
% arguments:  datarun - datarun struct (no required fields)
%            stimulus - struct specifying a stimulus to compare against datarun
%              source - string specifying the source of the provided stimulus
%
% outputs:    datarun - result of computation
%
%
% examples:
%
%  datarun = sync_stimulus(datarun, stimulus_from_stas,  'STAs');
%  datarun = sync_stimulus(datarun, stimulus_from_index, 'Index file');
%
%  2008  gauthier
%



% if there is no stimulus field in datarun, set datarun.stimulus to match the provided stimulus
if ~isfield(datarun,'stimulus')
    datarun.stimulus = stimulus;
    
else % if there is a datset.stimulus field, compare it to the provided stimulus

    % by default, assume all fields match.  this flag will change if a mismatch is found.
    mismatch = false;
    
    % initialize a struct that will be the union of datarun.stimulus and the provided stimulus
    union_stimulus = datarun.stimulus;

    % go through each of the provided fields
    provided_fields = fieldnames(stimulus);
    for ff = 1:length(provided_fields)
        
        % if the field doesn't exist in datarun.stimulus, add it
        if ~isfield(datarun.stimulus,provided_fields{ff})
            union_stimulus.(provided_fields{ff}) = stimulus.(provided_fields{ff});
            
        else % if it does exist, check for mismatch
            if ~isequal(stimulus.(provided_fields{ff}), datarun.stimulus.(provided_fields{ff}))
                %~all(stimulus.(provided_fields{ff}) == datarun.stimulus.(provided_fields{ff}))
                mismatch = true;
                break
            end
        end
    end

    % if there was a mismatch, print both versions of the stimulus and give an error
    if mismatch 
        fprintf('\n\ndatarun.stimulus\n\n')
        disp(orderfields(datarun.stimulus))
        fprintf('\nstimulus from %s\n\n',source)
        disp(orderfields(stimulus))
        warning('datarun.stimulus does not match stimulus from %s (see above).',source)
    else
        % otherwise, set datarun.stimulus to equal the union
        datarun.stimulus = union_stimulus;
    end
end

