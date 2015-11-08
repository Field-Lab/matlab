function r2 = compute_r2(prediction, data)

min_size = min(numel(data), numel(prediction));
prediction = prediction(1:min_size);
data = data(1:min_size);

if size(prediction,1)~=size(data,1)
    data = data';
end


sse = sum((prediction-data).^2); % variance between model and data
sst = sum((data-mean(data)).^2); % variance in the data
r2 = 1 - sse/sst;