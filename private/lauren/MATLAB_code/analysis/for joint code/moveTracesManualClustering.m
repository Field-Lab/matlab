function Data = moveTracesManualClustering(hObject, H, Data)

% function describing some callbacks of the gui manualCluster2D.m

% written by Lauren Hruby, SNL-E
% 2008-08-24

Data.prevDiscludedIndeces = Data.discludedIndeces;
Data.prevIncludedIndeces = Data.includedIndeces;

%getting information about the callback-initiating object and state of gui
if hObject == H.toLeftButton
    direction = 'left';
    if isempty(Data.discludedIndeces)
        errorh = errordlg('No traces to move left!', '', 'modal');
        uiwait(errorh);
        return
    end
elseif hObject == H.toRightButton
    direction = 'right';
    if isempty(Data.includedIndeces)
        errorh = errordlg('No traces to move right!', '', 'modal');
        uiwait(errorh);
        return
    end
else % swap button was pushed
    Data.includedIndeces  = Data.prevDiscludedIndeces;
    Data.discludedIndeces = Data.prevIncludedIndeces;
    Data.includedTraces  = Data.traces(Data.includedIndeces,:);
    Data.discludedTraces = Data.traces(Data.discludedIndeces,:);
    return
end


clusterType = get(H.selectionType, 'SelectedObject');

% clustering by PCA
if clusterType == H.usePCAButton
	if strcmp(direction,'left')
		[PCACoef, PCAScore] = princomp(Data.discludedTraces);
	else
		[PCACoef, PCAScore] = princomp(Data.includedTraces);
	end
	[PCx PCy] = PCChooser(PCAScore);

	%displays data and chosen PC plot for hand-clustering
	figure; set(gcf,'position',[50,50,800,1200]); title('lasso the points you want to move');
	[selx, sely, indeces] = lasso(PCAScore(:,PCx), PCAScore(:,PCy));
	close; close

% clustering by selecting traces directly
elseif clusterType == H.useTracesButton
	figh = figure; set(figh,'position',[50,50,800,1200]); title('lasso the traces you want to move')
	if strcmp(direction,'left')
        Data.discludedTraces
        size(Data.discludedTraces)
		indeces = lassoTraces(Data.discludedTraces);
	else
		indeces = lassoTraces(Data.includedTraces);
	end
	close

else
	errorh = errordlg('Selection type not specified!', '', 'modal');
	uiwait(errorh);
	return
end

% creates copy of current indeces so that removing elements doesn't mess up mapping from
% indeces to discluded/includedIndeces
includedIndecesCopy = Data.includedIndeces;
discludedIndecesCopy = Data.discludedIndeces;

if strcmp(direction,'left')
    for i = 1:length(indeces)
        moveToIncluded = discludedIndecesCopy(indeces(i));
        Data.discludedIndeces(Data.discludedIndeces == moveToIncluded) = [];
        Data.includedIndeces = [Data.includedIndeces moveToIncluded];
    end
else %if moving traces right
    for i = 1:length(indeces);
        moveToDiscluded = includedIndecesCopy(indeces(i));
        Data.includedIndeces(Data.includedIndeces == moveToDiscluded) = [];
        Data.discludedIndeces = [Data.discludedIndeces moveToDiscluded];
    end
end


Data.includedIndeces = sort(Data.includedIndeces);
Data.discludedIndeces = sort(Data.discludedIndeces);

%updates Data.included/discludedTraces
Data.includedTraces = Data.traces(Data.includedIndeces,:);
Data.discludedTraces = Data.traces(Data.discludedIndeces,:);

end