function poly = nyc2poly(nyc)
poly = cellfun(@segs2poly, nyc, 'UniformOutput', false);