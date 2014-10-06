
% assumes variable "results" exists

% print class names
for cc=1:length(results);fprintf('\n%d: %s',cc,results{cc}.name),end

class_list = [1 3 4 5 6 7 8 9 10]; %[1 2 3 5 6];%[2 7 9 10 12];

num_classes = length(class_list);
num_params = size(results{1}.cells,2);


% get number of unique cells for each params set
num_unique = zeros(num_classes,num_params+1);

for cc =1:num_classes
    % found by each params set
    num_unique(cc,1:num_params) = sum(results{class_list(cc)}.cells>0,1);
    % total unique
    num_unique(cc,num_params+1) = size(results{class_list(cc)}.cells,1);
end

% plot
figure(23);clf;barh(fliplr(num_unique'),'stacked')
if CELL
    labels = fo;
    labels{end+1} = 'union';
else
    for ll =1:size(fo,1);labels{ll} = fo(ll,:);end;labels{ll+1} = 'union';
end
set(gca,'YTick',1:(num_params+1),'YTickLabel',labels)
