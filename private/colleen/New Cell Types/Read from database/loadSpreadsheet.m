% extract info from database
clear

[~,TXT,~] = xlsread('/Volumes/Lab/Users/crhoades/Database spreadsheet 2015-12-1.xlsx', 'Dataruns');
[~,text_animal,~] = xlsread('/Volumes/Lab/Users/crhoades/Database spreadsheet 2015-12-1.xlsx', 'Piece');


gratings = strcmp(TXT, 'drifting-sinusoid');

[x,y] = find(gratings == 1);

non_empty_dates = find(~strcmp(TXT(:,1) ,''));
right_date = cell(length(x),2);
counter = 1;
for i = 1:length(x)
    all_dates_before = non_empty_dates(non_empty_dates <= x(i));
    right_date(i,1)= TXT(all_dates_before(end),1);
    right_date(i,2) = TXT(x(i),4);
    [a] = find(strcmp(text_animal(:,1), right_date(i,1))==1);
    if ~isempty(a)
        monkey_dates(counter,1) = TXT(all_dates_before(end),1);
        monkey_dates(counter,2) = TXT(x(i),4);
        counter = counter+1;
    end
    
end


