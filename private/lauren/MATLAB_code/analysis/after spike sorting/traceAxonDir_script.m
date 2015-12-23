
clear all

axonDirInfo = dataset_list_axon_directions();

for ii = 7:length(axonDirInfo)
    axonDirInfo(ii).direction = [];
    axonDirInfo(ii).endPts = {};
    axonDirInfo(ii).elecCoords = [];
    
    eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(axonDirInfo(ii).eiPath);
    
    for jj = 1:length(axonDirInfo(ii).axonTraceCells)
        ei = eiFile.getImage(axonDirInfo(ii).axonTraceCells(jj));
        ei = reshape(ei(1, 2:end, :), 64, []);
        
        [direction endPts elecCoords] = traceAxonDir(ei);
        axonDirInfo(ii).direction(end+1) = direction;
        axonDirInfo(ii).endPts{end+1} = endPts;
        if jj == 1;
            axonDirInfo(ii).elecCoords = elecCoords;
        end
        clear ei direction endPts elecCoords
    end
    eiFile.close()
    
    axonDirData = axonDirInfo(ii);
    savePath = axonDirInfo(ii).eiPath;
    fsInd = strfind(savePath, filesep);
    savePath = savePath(1:fsInd(end));
    
    %keyboard
    save([savePath 'axonDirData.mat'], 'axonDirData')
end

 %edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(elecResp.names.rrs_ei_path);