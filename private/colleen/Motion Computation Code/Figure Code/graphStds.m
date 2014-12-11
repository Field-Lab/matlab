%graphs for presentation:
close all
requests = {'ALL', 'ALLpooled'};

for j = 1:size(requests, 2)
    request = requests{j};
    % Options: BBMR, DBMR, BBML, DBML, ALL, tau, BBMRpooled
    if strcmp('BBMR', request)
        tag = 'Bright Bar Moving Right';
    elseif strcmp('DBMR', request)
        tag = 'Dark Bar Moving Right';
    elseif strcmp('BBML', request)
        tag = 'Bright Bar Moving Left';
    elseif strcmp('DBML', request)
        tag = 'Dark Bar Moving Left';
    elseif strcmp('ALL', request);
        tag = 'All Stimuli';
    elseif strcmp('tau', request);
        tag = 'Bright Bar Moving Right';
    elseif strcmp('BBMRpooled', request);
        tag = {'Bright Bar Moving Right'; 'ON Parasol pooled with ON Midget'};
    elseif strcmp('BBMLpooled', request);
        tag = {'Bright Bar Moving Left'; 'ON Parasol pooled with ON Midget'};
    elseif strcmp('DBMRpooled', request);
        tag = {'Dark Bar Moving Right'; 'ON Parasol pooled with ON Midget'};
    elseif strcmp('DBMLpooled', request);
        tag = {'Dark Bar Moving Left'; 'ON Parasol pooled with ON Midget'};
    elseif strcmp('ALLpooled', request);
        tag = {'All Stimuli'; 'ON Parasol pooled with ON Midget'};
    end
    
    
    data  = getData(request);
    if strcmp('BBMRpooled', request) || strcmp('BBMLpooled', request)|| strcmp('DBMRpooled', request)|| strcmp('DBMLpooled', request)|| strcmp('ALLpooled', request);
        meanOnPOnM = data(:,6);
        stdOnPOnM = data(:,7);
    else
        
        meanOnP(:,j) = data(:,6);
        stdOnP(:,j) = data(:,7);
        
        meanOffP(:,j) = data(:,8);
        stdOffP(:,j) = data(:,9);
        
        meanOnM(:,j) = data(:,10);
        stdOnM(:,j) = data(:,11);
        
        meanOffM(:,j) = data(:,12);
        stdOffM(:,j) = data(:,13);
        date = zeros(14,size(requests, 2));

    end
    
    filter(:,j) = data(:,5);
    delta(:,j) = data(:,4)*12;
    
    if strcmp('tau', request)
        tau(:,j) = data(:,14);
    end
    
    %  date = repmat('2007-08-24-4',8,1);
    for i = 1:sum(data(:,1) == 1976) %3/27/07
        ind = find(data(:,1) == 1976);
        %  date(ind(i),:) = '2007-03-27-1';
        date(ind(i),j) = 1;
        
    end
    
    
    %Normalize
    if strcmp('BBMRpooled', request) || strcmp('BBMLpooled', request) || strcmp('DBMRpooled', request)|| strcmp('DBMLpooled', request)|| strcmp('ALLpooled', request);
        stdOnPOnM = stdOnPOnM./delta(:,j);
        meanOnPOnM = convertMeanToDegSec(meanOnPOnM);
        
    else
        
        stdOnP(:,j) = stdOnP(:,j)./delta(:,j);
        stdOffP(:,j) = stdOffP(:,j)./delta(:,j);
        stdOnM(:,j) = stdOnM(:,j)./delta(:,j);
        stdOffM(:,j) = stdOffM(:,j)./delta(:,j);
        
        meanOnP(:,j) = convertMeanToDegSec(meanOnP(:,j));
        meanOffP(:,j) =  convertMeanToDegSec(meanOffP(:,j));
        meanOnM(:,j) =  convertMeanToDegSec(meanOnM(:,j));
        meanOffM(:,j) =  convertMeanToDegSec(meanOffM(:,j));
        
    end
    if ~(strcmp('BBMRpooled', request) || strcmp('BBMLpooled', request) || strcmp('DBMRpooled', request)|| strcmp('DBMLpooled', request)|| strcmp('ALLpooled', request))
        graphSpeed(meanOnP(:,j), meanOffP(:,j), meanOnM(:,j), meanOffM(:,j),stdOnP(:,j), stdOffP(:,j), stdOnM(:,j), stdOffM(:,j), delta(:,j), tag)
        graphByDate(stdOnP(:,j), stdOffP(:,j), stdOnM(:,j), stdOffM(:,j), date(:,j), tag)
        graphStimWidths(stdOnP(:,j), stdOffP(:,j), stdOnM(:,j), stdOffM(:,j), filter(:,j), tag)
        graphWithinRunVar(stdOnP(:,j), stdOffP(:,j), stdOnM(:,j), stdOffM(:,j), tag)
        if strcmp('tau', request)
            graphTau(stdOnP(:,j), stdOffP(:,j), stdOnM(:,j), stdOffM(:,j), tau(:,j), tag)
        end
    end
    
end

graphPooled(meanOnP, meanOffP, meanOnM, meanOffM,stdOnP, stdOffP, stdOnM, stdOffM,meanOnPOnM, stdOnPOnM, delta, tag);