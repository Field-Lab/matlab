% for each specified array id, plot the electrode positions


% choose array ids
array_ids = [100 500 1500 3500];


for array_id = array_ids

    % get java object
    electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(array_id);

    % note number of electrodes (skip electrode 0)
    num_electrodes = electrodeMap.getNumberOfElectrodes - 1;

    % get location of each electrode
    positions = zeros(num_electrodes,2);
    for pp = 1:num_electrodes
        positions(pp,1) = electrodeMap.getXPosition(pp);
        positions(pp,2) = electrodeMap.getYPosition(pp);
    end
   
    % plot electrodes
    figure(array_id);clf;
    plot(positions(:,1),positions(:,2),'.')
    axis equal; hold on
    
    % plot electrode IDs
    for ee=1:size(positions,1)
        text(positions(ee,1),positions(ee,2),num2str(ee),'Color',[.5 .5 0],'FontSize',12,...
            'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
    
end

