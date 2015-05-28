spikes_dir   = '/Volumes/Analysis/nora/NSEM/BlockedSpikes';

exp=4;
files1 = dir([spikes_dir '/2013-08-19-6/NSEM_mapPRJ']);
files2 = dir([spikes_dir '/2013-08-19-6/WN_mapPRJ']);


clear cell_id_1 cell_id_2 cell_id

count=0;
for i=1:length(files1)
    name=files1(i).name;
    if length(name)>20
        if strcmp(name(17:19),'OFF') && strcmp(name(20:22),'Par');
            count=count+1;
            cell_id_1(count)=str2double(name(24:(end-4)));
        elseif strcmp(name(17:18),'ON') && strcmp(name(19:21),'Par');
            count=count+1;
            cell_id_1(count)=str2double(name(23:(end-4)));
        end
    end
end

count=0;
for i=1:length(files2)
    name=files2(i).name;
    if length(name)>20
        if strcmp(name(17:19),'OFF') && strcmp(name(20:22),'Par');
            count=count+1;
            cell_id_2(count)=str2double(name(24:(end-4)));
        elseif strcmp(name(17:18),'ON') && strcmp(name(19:21),'Par');
            count=count+1;
            cell_id_2(count)=str2double(name(23:(end-4)));
        end
    end
end

count=0;
for i=1:length(cell_id_1)
   if ~isempty(cell_id_2(cell_id_2 == cell_id_1(i)))
      count=count+1;
      cell_id{count}=cell_id_1(i);
   end
end

all_cells{exp}=cell_id;

