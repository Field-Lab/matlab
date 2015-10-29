function mfr = max_firing_rate(raster, bin, duration)

x = bin/2:bin:duration-bin/2;
idx = find(~cellfun(@isempty,raster),1); % get the idx of 1st non-empty cell in raster
mfr = zeros(length(raster), length(raster{idx}));
for cc = 1:length(raster)
    if ~isempty(raster{cc})
        for tp = 1:length(raster{idx})
            h = hist(raster{cc}{tp}, x);
            mfr(cc, tp) = max(h);
        end
    end
end