

% notebooks directory
notes=dir('/Volumes/Data/Notebooks/20*');
expression='(\d+)-(\d+)-(\d+)';
voronoi_dates = cell(20,1);
cnt=1;
for i=1:length(notes)
    name=notes(i).name;
    mydate = regexpi(name,expression,'match');
    if ~isempty(mydate)
        
        myNotebook=['/Volumes/Data/Notebooks/',name];
        if myNotebook(end)=='d'
            myNotebook=[myNotebook,'/TXT.rtf'];
        end
        % attempt to read the notebook
        a=textread(myNotebook,'%s');
        
        for k=1:length(a)
            if ~isempty(regexpi(a{k},'voronoi')) %|| ~isempty(regexpi(a{k},'region'))
                voronoi_dates{cnt} = mydate;
                cnt = cnt+1;
                break;
            end
        end
    end
end
