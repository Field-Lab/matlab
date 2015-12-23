function bcf = combine_bcfs(bcfs)

% Check consistency as possible
for i = 1:(length(bcfs) - 1)
    if any(bcfs{i}.rois_x ~= bcfs{i+1}.rois_x)
        error('BCFs do not appear to have consistent ROIs');
    end
end


% Initialize with first
bcf = struct();
bcf.rois_x          = bcfs{1}.rois_x;
bcf.rois_y          = bcfs{1}.rois_y;
bcf.dll             = bcfs{1}.dll;
bcf.all_added_cones = bcfs{1}.all_added_cones;


% Add in the rest
for i = 2:length(bcfs)
    bcf.all_added_cones = [bcf.all_added_cones; bcfs{i}.all_added_cones];
    tocopy = bcf.dll == 0 & bcfs{i}.dll ~= 0;
    bcf.dll(tocopy) = bcfs{i}.dll(tocopy);
end