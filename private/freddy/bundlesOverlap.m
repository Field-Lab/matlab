function overlap = bundlesOverlap(shortBundle, longBundle)

overlapCount = 0;

bundleLen = size(shortBundle, 1);

for i = 1:bundleLen
    if ismember(shortBundle(i), longBundle)
        overlapCount = overlapCount + 1;
    end
end

if overlapCount/bundleLen > 0.8
    overlap = true;
else
    overlap = false;
end

end