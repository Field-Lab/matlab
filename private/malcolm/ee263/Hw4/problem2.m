lambda = [0.2 2.0 20];
ytilde = [y ; zeros(size(Edges,1), 1)];
yhat = nan(n, length(lambda));
A = eye(n);
F = nan(length(Edges), n);
for i = 1:length(Edges)
    F(i,:) = A(Edges(i,1),:) - A(Edges(i,2),:);
end

for i = 1:length(lambda)
    Atilde = [A ; sqrt(lambda(i)) * F];
    y_hat(:,i) = Atilde\ytilde; %#ok<SAGROW>
end

%laplacian_smoothing_plot(y_true,'True data, y_{true}',y,'Noisy data, y',0,'')
laplacian_smoothing_plot(y_hat(:,1),'y_{hat}, lambda = 0.2',y_hat(:,2),'y_{hat}, lambda = 2',y_hat(:,3),'y_{hat}, lambda = 20')