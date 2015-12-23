cd /snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002

fileIndex = 0;
files = dir;

%removes non-elecResp files/folders from list
for i = 1:length(files)
    i
    if ~isempty(strfind(files(i).name, 'elecResp'))
        fileIndex = fileIndex + 1;
        fileNames{fileIndex} = files(i).name; %#ok<SAGROW>
    end
end

errordlg('This script needs to be updated to reflect change in erfFitter to using maximum likelihood')



linErfParams = cell(length(fileNames), 1);
linErfErr = zeros(length(fileNames), 1);
linErfThresh = cell(length(fileNames), 1);

weibullParams = cell(length(fileNames), 1);
weibullError  = zeros(length(fileNames), 1);
weibullThresh = cell(length(fileNames), 1);

logErfParams = cell(length(fileNames), 1);
logErfErr = zeros(length(fileNames), 1);
logErfThresh = cell(length(fileNames), 1);

weibullConstrained = struct('params', cell(length(fileNames), 25), 'errors', cell(length(fileNames), 25), 'thresh', cell(length(fileNames), 25));
weibullConstrainedParams = zeros(length(fileNames), 25, 3);
weibullConstrainedErrors = zeros(length(fileNames), 25);


index = 0;
        
for i = 1:length(fileNames)
    temp = load([fileNames{i}]);
    elecResp = temp.elecResp;

    if isfield(elecResp, 'stimInfo') %if elecResp is loaded up
        nMovies = length(elecResp.stimInfo.movieNos);
        
        
        %checks to make sure values in successRates are correct
        for j = 1:nMovies
            mNum = elecResp.stimInfo.movieNos(j);
            if elecResp.stimInfo.nPulses(j) == 0 %value hasn't been set yet
                dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo, mNum, 99999);
                elecResp.stimInfo.nPulses(j) = size(dataTraces, 1);
            end
            if ~isempty(elecResp.analysis.type{j})
                successRate = sum(elecResp.analysis.latencies{j}~=0)/elecResp.stimInfo.nPulses(j);
                if elecResp.analysis.successRates(j) ~= successRate
                    disp(['warning: success rate for movie ' num2str(mNum) ' had to be corrected'])
                    elecResp.analysis.successRates(j) = successRate;
                end
            end
        end

        data = zeros(2, nMovies);
        data(1,:) = elecResp.stimInfo.stimAmps;
        data(2,:) = elecResp.analysis.successRates;
        lockedAmps = elecResp.analysis.finalized;
        for j = length(elecResp.stimInfo.stimAmps): -1: 1
            if isempty(elecResp.analysis.type{j})
                data(:,j) = [];
                lockedAmps(j) = [];
            end
        end
    end


    if data(2, end) > 0.5
        index = index + 1;
        
        index
        
        if index > 13
            keyboard
        end
            
        % linear-based
        data(1,:) = abs(data(1,:));
        [linErfParams{index} x linErfErr(index)] = erfFitter(data, 2, -1, 'makePlot', 0);
        linErfThresh{index} = -linErfParams{index}(2)/linErfParams{index}(1);
        

        [weibullParams{index} x weibullError(index)] = weibullFitter(data, 5, 1, 0, 'makePlot', 0);
        weibullThresh{index} = weibullParams{index}(2)*((-log(0.5))^(1/weibullParams{index}(1)))+weibullParams{index}(3); %check

        % log-based
        logdata = data;
        logdata(1,:) = log(abs(data(1,:)));
        [logErfParams{index} x logErfErr(index)] = erfFitter(logdata, 2, 1, 'makePlot', 0, 'lockedAmps', lockedAmps);
        logErfThresh{index} = -logErfParams{index}(2)/logErfParams{index}(1);
        
        %test weibull errors with different constrained values of k (a; first parameter)
        for j = 1:25
            [weibullConstrained(index, j).params x weibullConstrained(index, j).errors] =  weibullFitter(data, j, 1, 0, 'makePlot', 0, 'setParam', 1);
            weibullConstrained(index, j).thresh = weibullConstrained(index, j).params(2)*((-log(0.5))^(1/weibullConstrained(index, j).params(1)))+weibullConstrained(index, j).params(3); %check
            weibullConstrainedParams(index, j, :) = weibullConstrained(index, j).params;
            weibullConstrainedErrors(index, j) = weibullConstrained(index, j).errors;
        end
        
    end
end
finalIndex = index;

%plotting the results!!!

[sortedWeibullErrors sortIndeces] = sort(weibullError(1:35));

figure
plotColors = hsv(25);
hold on
for i = 1:25
    plot(weibullConstrainedErrors(sortIndeces, i), 'Color', plotColors(i, :));
end

plot(linErfErr(sortIndeces), 'k-', 'LineWidth', 2)
plot(weibullError(sortIndeces), 'r-', 'LineWidth', 2)
plot(logErfErr(sortIndeces), 'b-', 'LineWidth', 2)

hold off


