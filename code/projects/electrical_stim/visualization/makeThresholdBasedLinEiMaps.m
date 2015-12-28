function makeThresholdBasedLinEiMaps(dataPath, thrVec, varargin)
%function makeThresholdBasedLinEiMaps(dataPath, thrVec, [savePath], [soma or axon filter], [soma threshold filter if axon filter is specfied])
%Generate preparation surveys based on neuron main recording electrode thresholds, specified in second argument. 
%First argument contains Grind analysis output of visual stim
%Optional third argument is savePath, which gives base name to save image at (threshold value will be added)
%*you can put in a blank string to turn this off
%Optional fourth argument specifies whether to filter based on axonal threshold or soma threshold, or both.
%0 (default) means filter based on soma, 1 means filter based on axon. 2 sets a filter on the soma threshold
%Example: makeThresholdBasedEiMaps('/Volumes/Analysis/2015-11-09-10/data001/', [-20 -35 -50])
%Example: makeThresholdBasedEiMaps('/Volumes/Analysis/2015-11-09-10/data001/', [-20 -35 -50], '/Volumes/Analysis/2015-11-09-10/picture',1)
%Example: makeThresholdBasedEiMaps('/Volumes/Analysis/2015-11-09-10/data001/', [-20 -35 -50], '',1)
%Example: makeThresholdBasedEiMaps('/Volumes/Analysis/2015-11-09-10/data001/', [-20 -35 -50], '',2,-50)
%December 2015, Sasi Madugula
%*This function code is very messy because I've copied from surveyPrepartion.m and simply commented out the irrelevant bits

%parse varargin
if length(varargin) == 1 
    savStr = varargin{1}; 
    savFlg = 1; 
    filterFlg = 0;
elseif length(varargin) == 2;
    if isempty(varargin{1}); savFlg = 0; else; savStr = varargin{1}; savFlg = 1; end
    filterFlg = varargin{2};
elseif length(varargin) == 3
    if isempty(varargin{1}); savFlg = 0; else; savStr = varargin{1}; savFlg = 1; end
    filterFlg = 2;
    somaThr = varargin{3};
else
    filterFlg = 0;
    savFlg = 0;
end


%Make datarun variable
datarun  = load_data(dataPath); %this isn't doing its job on Windows for some reason, so doing it manually - Sasi
dataPathSplit = strsplit(dataPath, filesep);
datarun.names.rrs_neurons_path = [dataPath filesep dataPathSplit{end} '.neurons'];
datarun.names.rrs_ei_path = [dataPath filesep dataPathSplit{end} '.ei'];
datarun  = load_neurons(datarun);
datarun  = load_ei(datarun, 'all');
%datarun.names.rrs_sta_path = 'E:\2015-11-09-10\data001\data001.sta';
%datarun.names.rrs_params_path = 'E:\2015-11-09-10\data001\data001.params';
%datarun  = load_sta(datarun, 'load_sta', 'all');
%datarun  = load_params(datarun);

%loop over thresholds
for thr = thrVec
    if filterFlg
	if filterFlg == 1
	    [cellIdsToCheck, cellIndices, electrodes] = getLargeAmpSpikesAxon(datarun, thr); 
	else
	    [cellIdsToCheck, cellIndices, electrodes] = getLargeAmpSpikesAxon(datarun, thr, somaThr); 
	end
    else
	[cellIdsToCheck, cellIndices, electrodes] = getLargeAmpSpikes(datarun, thr); 
    end

    % Determine firing rate stability of these cells
    %*Lauren had this info in Data table, not sure what to do with it but will generate anyways. Might want to filter based on it
    N = length(cellIdsToCheck); 
    T = datarun.duration; % length of run
    binsize     =  5; % seconds
    binranges   = binsize:binsize:T; 
    firing_rate = zeros(N,1);
    firing_rate_vector = zeros(N,length(binranges)); 
    for i=1:N
	firing_rate(i) = length(datarun.spikes{cellIndices(i)})/T;
	bincounts = histc(datarun.spikes{cellIndices(i)},binranges); 
	firing_rate_vector(i,:) = bincounts/binsize; 
    end
    devs                = std(firing_rate_vector,[],2);
    binranges   = binranges; 
    firingRates = firing_rate_vector; 
     
    %goodCol = num2cell(false(N,1));
    %data = num2cell(cat(2,cellIdsToCheck',electrodes,firing_rate,devs));
    % toCat = num2cell(goodCol); 
    % data = {(cellIdsToCheck'),(electrodes),(firing_rate),(devs),(goodCol)};
    % data = {num2cell(cellIdsToCheck'),num2cell(electrodes),num2cell(firing_rate),num2cell(devs),num2cell(goodCol)};
    %data = [data goodCol]; 
    %set(handles.uitable1,'Data',data);
    % Update handles structure
    %guidata(hObject, handles);

    %plot linear axon traces

    % --- Executes on button press in poly_axon_traces.
    %function lin_axon_traces_Callback(hObject, eventdata, handles)
    fh = figure; 
    [xc, yc] = getElectrodeCoords512();
    %table_data = get(handles.uitable1,'Data');
    %colors = lines(size(table_data,1));
    colors = lines(N);
    nearby_axons = zeros(1, 512);
    nearby_somas = zeros(1, 512);
    nearby_range = 1; %measured in number of electrode distances, can be fractional
    for n = 1:1:N
	cellID = cellIdsToCheck(n);
	cellIndex = get_cell_indices(datarun, cellID);
	ei = datarun.ei.eis{cellIndex}'; % squeeze(ei(1,2:end,:))';
	eiAmps = max(ei)-min(ei);
	thresh = 4;
	[~,col,~] = find(eiAmps > thresh);
	aa = round(eiAmps(col))';
	yy = xc(col)';
	xx = yc(col)';
	
	wy = []; wx = [];
	% dumb way
	for i = 1:1:length(aa)
	    tmp = repmat(yy(i),aa(i),1);
	    wy = [wy; tmp];
	    tmp = repmat(xx(i),aa(i),1);
	    wx = [wx; tmp];
	end
     
	% vector of 1-D look-up table "x" points
	XI = linspace(min(xx),max(xx),max(round(abs(min(xx)-max(xx))/120),2));
	
	% obtain vector of 1-D look-up table "y" points
	YI = lsq_lut_piecewise( wx, wy, XI );
	sortaa = sort(aa,1,'descend'); 
	largestAmps = sortaa(1:2); 
	[~,IA,~] = intersect(aa,largestAmps); 
	
	COMx = 1/sum(largestAmps) * sum(xx(IA).*aa(IA)); 
	COMy = 1/sum(largestAmps) * sum(yy(IA).*aa(IA)); 
	   
    %     figure(fh); plot(XI,YI,'*-','Color',colors(n,:)); 
    % %     hold on; scatter(yc(row),xc(row),eiAmps(row)*6,colors(n,:),'filled');   % Plot eis
    % %     hold on; scatter(xx(IA),yy(IA),aa(IA)*6,colors(n,:),'filled'); % largest signals
    %     hold on; scatter(COMx,COMy,6*mean(largestAmps), colors(n,:),'filled');
    %     hold on; plot(XI,YI,'*-'); 
	
	figure(fh); plot(YI,XI,'*-','Color',colors(n,:)); 
    %     hold on; scatter(yc(row),xc(row),eiAmps(row)*6,colors(n,:),'filled');   % Plot eis
	close = zeros(512, 1);
	for ind = 1:(size(YI, 1)-1)
	    p1 = [YI(ind) XI(ind)]; p2 = [YI(ind+1) XI(ind+1)];
	    vec = [-(p1(2)-p2(2)) p1(1)-p2(1)];
	    vec = vec/norm(vec)*nearby_range*60;
	    c1 = p1+vec;
	    c2 = p1-vec;
	    c3 = p2-vec;
	    c4 = p2+vec;
	    xverts = [c1(1) c2(1) c3(1) c4(1)]; yverts = [c1(2) c2(2) c3(2) c4(2)];
	    contained = find(inpolygon(xc, yc, xverts, yverts));
	    
	    for elec = contained
		if l2pdist([p1; p2],[xc(elec), yc(elec)]) < 60
		    close(find(close == 0, 1, 'first')) = elec;
		end
	    end
	end
	close(close == 0) = [];
	close = unique(close);
	nearby_axons(close) = nearby_axons(close) + 1;
	close = zeros(512, 1);
	for ind = 1:size(xc, 2)
	    if pdist([COMy COMx; xc(ind) yc(ind)]) < 60
		close(find(close == 0, 1, 'first')) = ind;
	    end
	end
	close(close == 0) = [];
	close = unique(close);
	nearby_somas(close) = nearby_somas(close) + 1;
	
     
    %     hold on; scatter(xx(IA),yy(IA),aa(IA)*6,colors(n,:),'filled'); % largest signals
	hold on; scatter(COMy,COMx,6*mean(largestAmps), colors(n,:),'filled');
	text(double(COMy),double(COMx),num2str(cellID)); 
	
    %     figure; 
    %     hold on; scatter(xc(col),yc(col),eiAmps(col)*6,colors(n,:),'filled');   % Plot eis
    %     plot(YI,XI,'*-','Color',0.5*colors(n,:)); 
    %     hold on; scatter(COMy,COMx,6*mean(largestAmps), 0.5 * colors(n,:),'filled');
    %     text(double(COMy),double(COMx),num2str(cellID)); 
    end
    axis image; axis off; 
    if savFlg %if save option was specified --Sasi
	saveas(gcf, [savStr '_' num2str(abs(thr))], 'jpg')
    end
    %hold on; scatter(xc,yc,5,'black','filled'); 
     
    %{
    figure; scatter(xc,yc,300,nearby_axons,'filled'); colorbar; title('Axons');
    xlabel(['# axons within ' num2str(nearby_range) ' elec distance(s)']);
    axis image; axis off;
    set(findall(gca, 'type', 'text'), 'visible', 'on');
    figure; scatter(xc,yc,300,nearby_somas,'filled'); colorbar; title('Somas');
    xlabel(['# somas within ' num2str(nearby_range) ' elec distance(s)']);
    axis image; axis off;
    set(findall(gca, 'type', 'text'), 'visible', 'on');
    figure; 
    for x = 1:512
	text(xc(x)+20,yc(x)+20,num2str(x),'HorizontalAlignment','center', 'Color', 'white');
    end
    hold on;
    scatter(xc,yc,300,nearby_somas+nearby_axons,'filled'); colorbar; title('Both');
    xlabel(['# axons+somas within ' num2str(nearby_range) ' elec distance(s)']);
    axis image; axis off;
    set(findall(gca, 'type', 'text'), 'visible', 'on');
    %}
  %  handles.lin_nearby_axons = nearby_axons;
  %  handles.lin_nearby_somas = nearby_somas;
  %  guidata(hObject, handles);
 
end
