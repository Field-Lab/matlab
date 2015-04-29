function thresh = estimateBundleThresh(path, patternNo, visualThresh)

bundleMeans = getBundleVoltagesAStar(path, patternNo, false);
x = abs(bundleMeans(:, 2, 1));
y = abs(bundleMeans(:, 1, 1));
movieNos = bundleMeans(:, 3, 1);

index = find(movieNos == visualThresh, 1, 'first');
range = index:(index+5);
xToFit = x(range);
yToFit = y(range);

firstFew = y(1:10);
baseline = mean(firstFew);

coeffs = polyfit(xToFit, yToFit, 1);

fittedX = linspace(min(x), max(x), 200);
fittedY = polyval(coeffs, fittedX);

figure; scatter(abs(bundleMeans(:, 2, 1)), abs(bundleMeans(:, 1, 1)));
hold on; plot(fittedX, fittedY, 'r-', 'LineWidth', 3);

x_int = (baseline-coeffs(2))/coeffs(1);

thresh = x_int;

end