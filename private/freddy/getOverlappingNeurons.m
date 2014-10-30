function overlapping = getOverlappingNeurons(datarun, neurons, others)

to_do = others;

overlapping = zeros(size(neurons,2), size(others,2));

for n = 1:size(neurons,2)
    to_do(to_do == neurons(n)) = [];
    [x1, y1] = getCellEllipse(datarun, neurons(n));
    for test = to_do
        [x2, y2] = getCellEllipse(datarun, test);
        X = InterX([x1; y1], [x2; y2]);    
        if ~isempty(X)
            overlapping(n, find(overlapping(n, :)==0,1,'first')) = test;
        end
    end
end

end