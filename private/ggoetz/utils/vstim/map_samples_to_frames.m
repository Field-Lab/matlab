function map_samples_to_frames(fo, trefresh, nsamples, outputpath)
%MAP_SAMPLES_TO_FRAMES creates a mapping between samples and frames
%
%  MAP_SAMPLES_TO_FRAMES(FO, T_REFRESH, NSAMPLES, OUTPUTPATH) will
%  associate the samples 1:NSAMPLES of a recording with image indices. The
%  images are displayed in the order specified in FO, and appear at the
%  times defined in T_REFRESH (in sample units). 
%  For example, image k will have index FO(k) and will have appeared at
%  time T_REFRESH(k) in the recording. 
%  The resulting map is stored in binary format in the folder OUTPUTPATH.
%
%  Note: frame 0 is the special empty frame. All samples happening before
%  the first refresh time will be associated with frame 0.

if length(fo) ~= length(trefresh)
    error('Number of frames shown does not match refresh times.')
end

% Figure out the file name 
if exist(outputpath) %#ok<EXIST>
    assert(isdir(outputpath))
else
    mkdir(outputpath)
end
filename = split(outputpath, filesep);
filename = filename{end};
filename = sprintf('%s.stf', filename);

% Create progress bar
fprintf('Mapping samples to frames\n');
fprintf([repmat('.',1,80) '\n']);
ndots = 0;

% Create the mapping file
fid = fopen(fullfile(outputpath, filename), 'wb');
crefresh = 0;
for i = 1:nsamples
    % Progress bar
    if (i/nsamples)*80 > ndots
        fprintf('.');
        ndots = ndots + 1;
    end
    
    % Mapping
    if (crefresh < length(trefresh)) && (i >= trefresh(crefresh + 1))
        crefresh = crefresh + 1;
    end
    if crefresh
        fwrite(fid, uint32(fo(crefresh)), 'uint32');
    else
        fwrite(fid, uint32(0), 'uint32');
    end
end
fclose(fid);
fprintf('\nMapping done.\n');

end % map_samples_to_frames