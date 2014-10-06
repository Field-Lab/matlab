function plotBundleVoltage(patternNos, display)
bundleMeans = getBundleVoltagesAStar(patternNos, display);%getBundleVoltages(patternNos, display);
temp = load('~/Development/matlab-standard/private/lauren/MATLAB_code/analysis/dataset-specific/axonBundleThresholds_byPattern_2012_09_24_3_data008.mat');
axonBundleThresholds = temp.axonBundleThresholds_byPattern_2012_09_24_3_data008;
figure; xlabel('Stimulation amplitude (uA)'); ylabel('Average bundle voltage (mV)'); whitebg('black'); set(gcf, 'InvertHardCopy', 'off');
for patternIndex = 1:size(patternNos, 2)
    hold on; scatter(abs(bundleMeans(:, 2, patternIndex)), abs(bundleMeans(:, 1, patternIndex)));
    x=[axonBundleThresholds(patternNos(patternIndex)),axonBundleThresholds(patternNos(patternIndex))];
    y=[0,300];
    hold on; plot(x,y); hold off;
%     f = @(F,x) (1 +exp(-F(1)*(x - F(2)))).^(-1); % sigmoid
%     F_fitted = nlinfit(abs(bundleMeans(:, 2, patternIndex)),abs(bundleMeans(:, 1, patternIndex)),f,[1 1]);              
%     y = f(F_fitted,abs(bundleMeans(:, 2, patternIndex)));
%     hold on; plot(abs(bundleMeans(:, 2, patternIndex)), y); hold off;
end

x=[0,4.5];
y=[29.3796,29.3796];
hold on; plot(x,y); hold off;

x=[0,4.5];
y=[16.05148356,16.05148356];
hold on; plot(x,y); hold off;

singleAxonVoltage = 16.05148356;
ylims = ylim;
axonupper = ylims(2)/singleAxonVoltage;

legend(num2str(patternNos'), -1);

set(gca,'Box','off');   %# Turn off the box surrounding the whole axes
axesPosition = get(gca,'Position');          %# Get the current axes position
hNewAxes = axes('Position',axesPosition,'Color','none','YLim',[0 axonupper],'YAxisLocation','right','XTick',[],'Box','off');                %#   ... and no surrounding box
ylabel(hNewAxes,'Approx. # of axons');
end