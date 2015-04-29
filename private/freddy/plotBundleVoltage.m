function plotBundleVoltage(path, patternNos, display)
bundleMeans = getBundleVoltagesAStar(path, patternNos, display);%getBundleVoltages(patternNos, display);

figure; xlabel('Stimulation amplitude (uA)'); ylabel('Average bundle voltage (mV)'); set(gcf, 'InvertHardCopy', 'off');
for patternIndex = 1:size(patternNos, 2)
    hold on; scatter(abs(bundleMeans(:, 2, patternIndex)), abs(bundleMeans(:, 1, patternIndex)));
end

singleAxonVoltage = 16.05148356;

%x=[0,4.5];
%y=[singleAxonVoltage,singleAxonVoltage];
%hold on; plot(x,y); hold off;

ylims = ylim;
axonupper = ylims(2)/singleAxonVoltage;

legend(num2str(patternNos'), -1);

set(gca,'Box','off');   %# Turn off the box surrounding the whole axes
axesPosition = get(gca,'Position');          %# Get the current axes position
hNewAxes = axes('Position',axesPosition,'Color','none','YLim',[0 axonupper],'YAxisLocation','right','XTick',[],'Box','off');                %#   ... and no surrounding box
ylabel(hNewAxes,'Approx. # of axons');
end