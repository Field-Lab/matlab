function [raster_sum, n_idx] = get_ndirection_raster(raster, pdirection)

idx = find(~cellfun(@isempty,raster),1); % get the idx of 1st non-empty cell in raster
spn = size(raster{idx}, 1);
tpn = size(raster{idx}, 2);
dirn = size(raster{idx}, 3);
pdirection(pdirection<0) = pdirection(pdirection<0) + 2*pi;
p_idx = round(pdirection/(2*pi/dirn))+1;
p_idx(p_idx == dirn+1) = 1;
n_idx = p_idx+4;
n_idx(n_idx>dirn) = n_idx(n_idx>dirn)-dirn;

raster_sum = cell(length(raster), 1);
for cc = 1:length(raster)
    if ~isempty(raster{cc})
        for sp = 1:spn
            for tp = 1:tpn
                raster_sum{cc}{sp, tp} = sort(cell2mat(squeeze(raster{cc}(sp, tp, n_idx(cc), :))));
            end
        end
    end
end

end