function corr_coef = compute_corr(prediction, data)

min_size = min(numel(data), numel(prediction));
prediction = prediction(1:min_size);
data = data(1:min_size);

if size(prediction,1)<size(prediction,2)
    prediction = prediction';
end

if size(data,1)<size(data,2)
    data = data';
end

corr_coef = corr(prediction, data);