function [ BadElectrodes ] = getBadElectrodes( DataPath, PatternRange, MovieNumber)

    maxs = zeros(512,1);
    for PatternNumber = PatternRange
        [DataTraces0,~,~] = NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);
        s = reshape(mean(DataTraces0),512,140);
        tmp = mean(s(:,135:end),2);
        s = s - repmat(tmp,1,140);
        st = s(:,7:40);
        [minimus,~] = min(st,[],2);
        values = -minimus;

        maxs = max([values maxs],[],2);
    end
    bad = (maxs < 100).*((1:512)');
    BadElectrodes = intersect(bad(bad ~= 0),PatternRange);
end

