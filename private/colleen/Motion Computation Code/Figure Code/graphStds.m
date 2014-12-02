%graphs for presentation:
close all
request = 'DBML';


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
end

data  = getData(request);

meanOnP = data(:,6);
stdOnP = data(:,7);

meanOffP = data(:,8);
stdOffP = data(:,9);

meanOnM = data(:,10);
stdOnM = data(:,11);
 filter = data(:,5);
 delta = data(:,4)*12;
 date = zeros(14,1);
%  date = repmat('2007-08-24-4',8,1);
 for i = 1:sum(data(:,1) == 1976) %3/27/07
     ind = find(data(:,1) == 1976);
%  date(ind(i),:) = '2007-03-27-1';
  date(ind(i),:) = 1;

 end
 
 %Normalize
stdOnP = stdOnP./delta;
stdOffP = stdOffP./delta;
stdOnM = stdOnM./delta;



graphSpeed(stdOnP, stdOffP, stdOnM, delta, tag)
graphByDate(stdOnP, stdOffP, stdOnM, date, tag)
graphStimWidths(stdOnP, stdOffP, stdOnM, filter, tag)
graphWithinRunVar(stdOnP, stdOffP, stdOnM, tag)