[x, y] = getElectrodeCoords61();
x = x';
y = y';
xy = [x y];

adj_mat_61 = cell(64,1);
for n = 1:64
    v = xy(n, :);
    if ~isnan(v)
        for a = 1:64
            if ~isnan(a)
                e = xy(a, :);
                dist = pdist([v;e]);
                if abs(pdist([v;e]) - 2) < 1
                    adj_mat_61{n} = [adj_mat_61{n} a];
                end
            end    
        end  
    else
        adj_mat_61{n} = NaN;
    end    
end
