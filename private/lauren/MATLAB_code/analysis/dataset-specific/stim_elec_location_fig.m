
%parameters for all datasets
%clear all

%analysis parameters
params.recalcAll = false; %redo all response curve fits and standard deviation calculations
params.nBootstrapReps = 100;
params.maxSuccRateCutoff = 0.4;
params.removeOffArrayCells = true; params.cutoffMult = 0.5;
params.exclude30ArrayData = true;
params.includeDataset = [true true true true false false false];

params.cellTypes = {'onPar', 'offPar', 'onMidg', 'offMidg', 'sbc'};

cellInfo = load_analysis_all_paper1_cells(params);


% ei soma location parameters
sigThresh = 0.1; %minimum abs(signal) required for an electrode to be included in COM calculation


%colors!
blue = [50 70 247]/255; %blue
rust = [.8 .05 0.05];  %rust
grass = [90 156 0]/255; %pale grass
salmon = [255 124 59]/255;  %salmon

colors(1,:) = 1-(1-blue)*0.7; %ON parasol
colors(2,:) = blue*0.7; %OFF parasol
colors(3,:) = 1-(1-rust)*0.7; %ON midget
colors(4,:) = rust*0.7; %OFF midget
colors(5,:) = grass; %SBC

%plotting details
[xCoords yCoords] = getElectrodeCoords61();

arrayWidth = max(xCoords)*2*1.2;
arrayHeight = max(yCoords)*2*1.2;

makePlot = true;
checkForAxons = false;


%%

allDists = [];

checkAxonFig = figure('position', [100 100 400 350]);
figure; plotAxes = axes();
for ii = 1:length(cellInfo)
    for jj = 1:length(cellInfo{ii})
        ei = cellInfo{ii}(jj).ei;
        
        % plot ei in array layout to look for axon signals
        
        %eiSigRange = [min(min(ei([1:8 10:24 26:56 58:64], :))) max(max(ei([1:8 10:24 26:56 58:64], :)))];
        %peakAmp = max(abs(eiSigRange));
        
        %calculate peak (absolute) signal on each electrode and normalize
        eiAmps = max(abs(ei),[],2);
        eiNegAmps = max(-1*(ei),[],2);

        peakElec = find(eiAmps == max(eiAmps));
        
        if peakElec ~= find(eiNegAmps == max(eiNegAmps));
            disp('WARNING: peak negative electrode and peak abs electrode are not the same!')
            keyboard
        end

        ei = ei/max(eiAmps);
        eiAmps = eiAmps/max(eiAmps);
                
        eiTmp = zeros(1,65,size(ei,2)); eiTmp(1,2:end,:) = ei;
        getCOM(eiTmp, eiAmps, sigThresh, makePlot, checkForAxons, cellInfo{ii}(jj).shortName, 'plotAxes', plotAxes,...
            'plotFigAxonCheck', checkAxonFig, 'axonElecs', cellInfo{ii}(jj).axonElecs);
        
        axes(plotAxes); hold on
        
        e1 = peakElec;
        e2 = cellInfo{ii}(jj).stimElec;
        
        plot(xCoords(e1), yCoords(e1), 'rx', 'markersize', 10)
        plot(xCoords(e2), yCoords(e2), 'b*')
        
        %pause
        
        xDist = abs(xCoords(e1)-xCoords(e2))*30*60/51.9615; %horizontal distance between electrode columns should be 51.9615 but are actually 60
        yDist = abs(yCoords(e1)-yCoords(e2))*30;
        
        cellInfo{ii}(jj).dist = norm([xDist yDist]);
        
        allDists = [allDists cellInfo{ii}(jj).dist];
                
        %clf(gcf)
        
%        figure('position', [100 100 400 350])
%        
%         for i = 1:64
%             if ~(i==9||i==25||i==57)
%                 axes('position', [xCoords(i)/arrayWidth + 0.45, yCoords(i)/arrayHeight + 0.45, 0.1, 0.1])
%                 plot(squeeze(ei(1,i+1,:)), 'k')
%                 set(gca, 'ylim', [-10 8], 'xlim', [1 81], 'xtick', [], 'ytick', [])
%             end
%         end
    end
end






%%



