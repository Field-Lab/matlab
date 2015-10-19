% spatial sensitivity plots comparing algorithm and human
cd('/Users/gomena/Research/EJBigData/Datasetsvisitjun15')
load Output216dif
	cd('/Users/gomena/Research/Git/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM')
% % cells ids from /Users/gomena/Research/EJBigData/Datasetsvisitjun15/2015-04-14-0/data000/data000.ei
% cellIds = [1 392 633 917 1308 1699 2345 3049 3466 3947 4655 6873];
[xc,yc] = getElectrodeCoords512(); 
codebase_path = matlab_code_path; 
cmap = load([codebase_path 'code/projects/electrical_stim/'...
    'resources/redtealcmap.mat']); 

% cells ids from /Users/gomena/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/data000/data000.ei
dataPath = '/Users/gomena/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/data000/data000'; 
cellIds = [3457 4058 4656 2796 5123 5748 6439];

% Analyze all electrodes that recorded voltages over some threshold. 
% Load EI. 
datarun  = load_data(dataPath);
datarun  = load_neurons(datarun);
datarun  = load_sta(datarun, 'load_sta', 'all');
datarun  = load_params(datarun);
datarun  = load_ei(datarun, cellIds,'keep_java_ei','false');

ei_thresh = 10; 
%%
figure(100); 
figure(300);
figure(200); scatter(xc,yc,3,'k','filled');    axis off; axis image;
line([0 0],[min(yc) max(yc)],'Color','k');
line([min(xc)/2 min(xc)/2],[min(yc) max(yc)],'Color','k');
line([max(xc)/2 max(xc)/2],[min(yc) max(yc)],'Color','k');


for n = 1:length(cellIds)
    thresh_quad1 = zeros(2,512);
    thresh_quad2 = zeros(2,512);
    thresh_quad3 = zeros(2,512);
    thresh_quad4 = zeros(2,512);
    neuronId = cellIds(n);
    ei = datarun.ei.eis{get_cell_indices(datarun,neuronId)};
    eiAmps = max(ei,[],2) - min(ei,[],2);
    if n == 1
        plotCoords = true;
    else
        plotCoords = false;
    end; 
    eiContour_wLinFit(eiAmps,'linFitThresh',5,'figureNum',100,'plotCoords',plotCoords);
    eiContour_wLinFit(eiAmps,'linFitThresh',5,'figureNum',300,'plotCoords',plotCoords);
    
    % Find ei amps greater than a particular threshold
    elecs = find(eiAmps > ei_thresh); 
    figure(200); 
    hold on; scatter(xc(elecs),yc(elecs),eiAmps(elecs),'r','filled'); 
 
    % Determine the quadrant for each electrode in a given cell .  
    
    if ~isempty(find(xc(elecs) < min(xc/2), 1))
        fpath = '/Users/gomena/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/data003/';
        [aa bb thresh_quad1 cc] = genActThreshSpatialMaps_Gonzalo(fpath, neuronId,Output216,...
            'dispCheckedElecs',1);
    end

    if ~isempty(find(xc(elecs) > min(xc/2) & xc(elecs) < 0, 1));
        fpath = '/Users/gomena/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/data004/';
        [aa bb thresh_quad2 cc] = genActThreshSpatialMaps_Gonzalo(fpath, neuronId,Output216,...
            'dispCheckedElecs',1);
    end
    
    if ~isempty(find(xc(elecs) < max(xc/2) & xc(elecs) > 0, 1))
        fpath = '/Users/gomena/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/data005/';
        [aa bb thresh_quad3 cc] = genActThreshSpatialMaps_Gonzalo(fpath, neuronId,Output216,...
            'dispCheckedElecs',1);

    end
    if ~isempty(find(xc(elecs) > max(xc/2), 1))
        fpath = '/Users/gomena/Research/EJBigData/Datasetsvisitjun15/2012-09-24-3/data006/';
        [aa bb thresh_quad4 cc] = genActThreshSpatialMaps_Gonzalo(fpath, neuronId,Output216,...
            'dispCheckedElecs',1);

    end
   cd('/Users/gomena/Research/GIT/ChichilniskyGIT/matlab/code/projects/electrical_stim/spike_sorting_GM/Manuscript/ResultsFigures/')
    thresholds = thresh_quad1 + thresh_quad2 + thresh_quad3 + thresh_quad4;
    titles={'Activation Thresholds (\muA), Human spike sorting','Activation Thresholds (\muA), Algorithm'};
    names={'HumanActSeveralNeurons','AlgActSeveralNeurons'};
    for i=1:2
    thresholds(i,thresholds(i,:)>6) = 0; 
    thresholds(i,thresholds(i,:)<0) = 0; 
    
    idx = find(thresholds(i,:));
    figure(100+200*(i-1)); 
    hold on; scatter(xc(idx),yc(idx),150, thresholds(i,idx),'filled');
   title(titles{i},'fontsize',1)
    colormap(cmap.m); caxis([0.5 4.5]);
    colorbar
    print(names{i},'-dpng')
     print(names{i},'-depsc2')
     
    end
end

