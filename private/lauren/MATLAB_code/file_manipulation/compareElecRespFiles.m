%% checks whether elecResp files in different directories are identical to eachother

clear all

% directories containing elecResp files
dir1 = '/marte/snle/lab/Experiments/Array/Analysis/2008-08-27-2/';
dir2 = '/marte/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/';

d1 = dir(dir1);
d2 = dir(dir2);

d2checked = false(length(d2),1);

for ii = 1:length(d1)
    if ~isempty(strfind(d1(ii).name, 'elecResp'))
        foundMatch = false;
        for jj = 1:length(d2)
            if strcmpi(d1(ii).name, d2(jj).name) %check if data in elecResp files is identical
                disp('checking...')
                elecResp1 = load([dir1 filesep d1(ii).name]);
                elecResp2 = load([dir2 filesep d2(jj).name]);
                keyboard
                if ~isequal(elecResp1, elecResp2)
                    disp([d1(ii).name ' does not contain identical data in both directories'])
                    if ~isequal(elecResp1.elecResp.analysis, elecResp2.elecResp.analysis)
                        disp([d1(ii).name ' does not contain identical analysis in both directories'])
                    end
                    keyboard
                end
                
                foundMatch = true;
                d2checked(jj) = true;
                
                clear elecResp1 elecResp2
                break
            end
        end
        if ~foundMatch
            disp([d1(ii).name ' exists only in the first directory'])
        end
    end
end

for jj = 1:length(d2)
    if ~isempty(strfind(d2(jj).name, 'elecResp')) && ~d2checked(jj)
        disp([d2(jj).name ' exists only in the second directory'])
    end
end