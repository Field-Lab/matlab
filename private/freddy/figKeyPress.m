%%
function figKeyPress(hObject, evt, handles)
cla;
act_max = 1;
act_min = 0;

handles = guidata(hObject);
datarun = handles.datarun;
neurons = handles.neurons;
activationThresh = handles.activationThresh;
responseCurves = handles.responseCurves;
currentIndex = handles.currentIndex;

if strcmpi(evt.Key, 'uparrow')
    currentIndex = currentIndex+1;
elseif strcmpi(evt.Key, 'downarrow')
    currentIndex = currentIndex-1;
end

handles.currentIndex = currentIndex;
guidata(hObject, handles);

responseProbs = responseCurves(:,:,currentIndex);

tmp = figure;
rgbVals = colormap(tmp,gray);
close(tmp);
for n = 1:1:size(handles.neurons,2)
    currentActThresh_single = activationThresh(n,:);
    currentProb = responseProbs(n,:);
    [~,J,val] = find(currentActThresh_single);
    minThresh = min(val);
    bestStimElec = J(find(val == min(val)));
    [~,~,val] = find(currentProb);
    if ~isempty(val)
        [rgbIndex, ind] = max(val);
        rgbIndex = round(rgbIndex * size(rgbVals,1));
        if rgbIndex<1; rgbIndex = 1; elseif rgbIndex>size(rgbVals,1); rgbIndex = size(rgbVals,1); end
        
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color',rgbVals(rgbIndex,:));
    else
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color','k');
    end
    
    %disp(['prob for neuron ' num2str(neurons(n)) ' = ' num2str(max(val)) ' at electrode ' num2str(bestStimElec) ' (single electrode stim)']);
end
axis image; axis off;
colorbar; colormap(rgbVals); caxis([act_min act_max]);
title('parasol response prob, single electrode stimulation');
end