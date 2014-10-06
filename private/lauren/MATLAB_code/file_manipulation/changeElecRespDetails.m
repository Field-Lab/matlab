clear all
%cd /Volumes/Palace/Analysis/Lauren/2008-08-26-0/data008
%neuronNo = 257;
%patternNos = 111:220;

%for i = patternNos

%cd /Analysis/lauren/2011-06-24-5/data015/
cd /snle/lab/Experiments/Array/Analysis/2011-08-04-5/data004/

fileIndex = 0;
files = dir;

%removes non-elecResp files/folders from list
for i = 1:length(files)
    i
    if ~isempty(strfind(files(i).name, 'elecResp')) && isempty(strfind(files(i).name, '._'))
        fileIndex = fileIndex + 1;
        fileNames{fileIndex} = files(i).name; %#ok<SAGROW>
    end
end

clear fileIndex files i

%%
for xx = 1:length(fileNames)
    
    xx
    disp(['examining ' fileNames{xx}])
    load(fileNames{xx})


    elecResp.names.data_path = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data004/';
    elecResp.names.artifact_path = '';
    elecResp.names.savePath = elecResp.names.data_path;
    elecResp.names.rrs_ei_path = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data000/data000.ei';
    elecResp.names.rrs_params_path = '/snle/lab/Experiments/Array/Analysis/2011-08-04-5/data000/data000.params';
    
%     elecResp.names.data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data008/';
%     elecResp.names.rrs_ei_path = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data007-NW/data007-NW.ei';
%     elecResp.names.rrs_params_path = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data007-NW/data007-NW.params';
%     elecResp.names.savePath = elecResp.names.data_path;
        
    %elecResp.names.rrs_params_path = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data003/data003.params';
    %rrs_short_name = '2010-11-22-1/data005_n286_p20_f10';
    %if ~strcmp(elecResp.names.rrs_params_path, '/Analysis/lauren/2008-08-27-2/data001-lh/data001-lh.params')
        

        %load(['elecResp_n' num2str(neuronNo) '_p' num2str(i) '.mat'])

        %elecResp.names.data_path = '/Volumes/Brokedown/Analysis/Lauren/2009-09-05-0/data009/';
        %elecResp.names.data_path = '/Analysis/lauren/2008-08-27-2/data002/';
        %elecResp.names.savePath = elecResp.names.data_path;
        %elecResp.names.data_path = 'C:\Documents and Settings\Lauren Hruby\desktop\lab\data014_015\';

        %elecResp.names.rrs_ei_path = '/Analysis/lauren/2008-08-27-2/data001-lh/data001-lh.ei';
        %elecResp.names.rrs_params_path = '/Analysis/lauren/2008-08-27-2/data001-lh/data001-lh.params';
        
        
        % to fix problem of saving extra variables into elecResp file
        %     currentVars = who;
        %     for j = 1:length(currentVars)
        %         if ~strcmpi(currentVars{j}, 'elecResp') && ~strcmpi(currentVars{j}, 'fileNames')...
        %                 && ~strcmpi(currentVars{j}, 'currentVars') && ~strcmpi(currentVars{j}, 'xx')
        %             clear(currentVars{j})
        %         end
        %     end


    %%%%%%%%%%  updates for elecResp files created by older versions of the
    %%%%%%%%%%  analysis software %%%%%%%%%%

    %     %%%%%%%%%%%%%%%% use if elecResp.cells.all doesn't exist %%%%%%%%%%%%%%%%%%%%%%%%
    %     electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
    %     channelsToUse = electrodeMap.getAdjacentsTo(elecResp.cells.recElec, 1)';
    % 
    %     for i = 1:length(elecResp.cells.goodElecs)
    %         if ~any(channelsToUse == elecResp.cells.goodElecs(i))
    %             channelsToUse = [channelsToUse elecResp.cells.goodElecs(i)]; %#ok<AGROW>
    %         end
    %     end
    %     elecResp.cells.all = findNeuronsAboveThresh(elecResp.names.rrs_ei_path,...
    %         elecResp.names.rrs_params_path, channelsToUse, 10);
    % 
    %     % adds main and active cells to 'all' list if not in it already, and checks to make sure each exists
    %     % in the params file
    %     for i = 1:length(elecResp.cells.active{1})
    %         if ~any(datarun.cell_ids == elecResp.cells.active{1}(i))
    %             error('Neuron ID specified doesn''t exist in the params file');
    %         elseif ~any(elecResp.cells.all == elecResp.cells.active{1}(i))
    %             elecResp.cells.all = [elecResp.cells.all elecResp.cells.active{1}(i)];
    %         end
    %     end
    %     if ~any(elecResp.cells.all == elecResp.cells.main)
    %         elecResp.cells.all = [elecResp.cells.all elecResp.cells.main];
    %     end

    
%     if ~isfield(elecResp.cells, 'mainEI')
% 
%         %%%%%%%%%%%%%%%% use if elecResp.cells.mainEI and .allEIs doesn't exist %%%%%%%%%%%%%%%%%%%%%%%%
%         %retrieves ei data for cells in elecResp.cells.all
%         eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(elecResp.names.rrs_ei_path);
%         elecResp.cells.mainEI = eiFile.getImage(elecResp.cells.main);
%         elecResp.cells.mainEI = reshape(elecResp.cells.mainEI(1, 2:end, :), 64, []);
%         
%         elecResp.cells.allEIs = cell(length(elecResp.cells.all),1);
%         for k = 1:length(elecResp.cells.all)
%             elecResp.cells.allEIs{k} = eiFile.getImage(elecResp.cells.all(k));
%             elecResp.cells.allEIs{k} = reshape(elecResp.cells.allEIs{k}(1, 2:end, :), 64, []);
%         end
%         clear eiFile
%         
%         
%         %%%%%%%%%%%%%%% use if .stimInfo.nPulses doesn't exist or is suspected to be incorrect %%%%%%%%%
%         elecResp.stimInfo.nPulses = zeros(1, length(elecResp.stimInfo.movieNos));
%         % determines number of pulses in each movie
%         for i = 1:length(elecResp.stimInfo.movieNos)
%             dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
%                 elecResp.stimInfo.movieNos(i), 99999);
%             elecResp.stimInfo.nPulses(i) = size(dataTraces, 1);
%         end
%         
%         %     %%%%%%%%%%%%%%% use if .cells.active and .analysis.otherLatencies don't exist%%%%%%%%%%%%%%%
%         %     %%%%%%%%%%%%%%% warning: be careful not to overwrite existing analysis!!! %%%%%%%%%%%%%%%%%%
%         
%         %     if ~isnumeric(elecResp.cells.active)
%         %
%         %     elseif ~all(elecResp.cells.active == elecResp.cells.main)
%         %
%         %     else
%         %         keyboard
%         %         elecResp.analysis.otherLatencies = cell(length(elecResp.stimInfo.movieNos), 1);
%         %         elecResp.cells.active = cell(length(elecResp.stimInfo.movieNos),1);
%         %     end
%         
%         %save(['elecResp_n' num2str(neuronNo) '_p' num2str(i) '.mat'], 'elecResp')
%         
%         %end
%         save([fileNames{xx}], 'elecResp')
%     else
%         keyboard
%     end

    save(fileNames{xx}, 'elecResp')
    
    
        
    clear elecResp
end
