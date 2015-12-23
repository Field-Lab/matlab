function s = form_cell_prefix(cellids)

s = '';

if (isempty(cellids))
    return;
end

for j=1:length(cellids)-1
    
    s = strcat(s,num2str(cellids(j)),'_');
    
end

s = strcat(s,num2str(cellids(end)));