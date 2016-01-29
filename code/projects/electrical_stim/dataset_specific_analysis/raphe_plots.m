function raphe_plots()
if 0
% Plot EIs from the Raphe recording showing EIs along the hair part.
dataPath = '/Volumes/Analysis/2015-10-29-5/data000/data000';%'/Volumes/Analysis/2012-09-24-3/data007/data007';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 'all');

% neuronIDs = [467 2794 3323 4216 4953 6181 6290 7638];
neuronIDs = [467 6181 2480 3138 4503 ];%932
badElecs = [];
plotContours(datarun,neuronIDs,badElecs,301);

%% Try a single cell
plotContours(datarun,932,badElecs,202);
title('neuron 932');
%% Show activation plots with electrical stimulation to show that we can analyze cells
plotArtifactSubtractedData('/Volumes/Analysis/2015-12-18-3/data001/',...
    'elecResp_n2867_p184.mat', [2.4 3.5] );
plotArtifactSubtractedData('/Volumes/Analysis/2015-12-18-3/data001/',...
    'elecResp_n2867_p192.mat', [0.8 1.3] );
plotArtifactSubtractedData('/Volumes/Analysis/2015-12-18-3/data001/',...
    'elecResp_n3106_p208.mat', [0.8 1.3] );
plotArtifactSubtractedData('/Volumes/Analysis/2015-12-18-3/data001/',...
    'elecResp_n5208_p348.mat', [0.4 0.8] );
end
%% Show a group of EIs from a typical experiment that all go one direction.
% dataPath = '/Volumes/Acquisition/Analysis/2015-12-18-1/data000/data000';%'/Volumes/Analysis/2012-09-24-3/data007/data007';
% datarun = load_data(dataPath);
% datarun = load_neurons(datarun);
% datarun = load_ei(datarun, 'all');
% neuronIDs = [3138 3902 4744  7038]; %5161 5585
dataPath = '/Volumes/Analysis/2015-10-29-1/streamed/data009/data009';
dataPath = '/Volumes/Analysis/2015-10-29-1/data009/data009';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_ei(datarun, 'all');
neuronIDs = [2071 3631 4292 5086 5807]; 

badElecs = [];
plotContours(datarun,neuronIDs,badElecs,400);

    function plotContours(datarun,neuronIDs,badElecs,figureNum)
        for n = 1:length(neuronIDs)
            neuronID = neuronIDs(n);
            EI = datarun.ei.eis{get_cell_indices(datarun, neuronID)};
            eiAmps = max(EI,[],2) - min(EI,[],2);
            for b = 1:length(badElecs)
                badElectrode = badElecs(b);
                cluster = getCluster519(badElectrode);
                newamps = mean(eiAmps(cluster(2:end)));
                eiAmps(badElectrode) = newamps;
            end
            figure(10); hold on; 
            
            if n==1
                showElectrodes = true;
                scatter(datarun.ei.position(:,1),datarun.ei.position(:,2),...
                    8,0.8*[1 1 1],'filled'); 
            else
                showElectrodes = false;
            end
            idx = find(eiAmps > 10); 
            scatter(datarun.ei.position(idx,1),datarun.ei.position(idx,2),eiAmps(idx),'filled'); 
            axis image; axis off; 
            eiContour_wLinFit(eiAmps,'threshForContours',4,'numContourLevels',8,'numElecs',519,...
                'figureNum',figureNum,'linFitThresh',4,'showSoma',0,'plotCoords',showElectrodes);
        end
    end

end