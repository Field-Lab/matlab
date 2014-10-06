
function paired_cells=pick_neighbor_cells(mean, parasol_struct, datarun_vision)

% find on and off parasol cell types

% for on and off parasols
 for j=1:2
     distance=zeros(parasol_struct.NumCells{j},1);
     for i=1:parasol_struct.NumCells{j}
         distance(i)=norm(datarun_vision{parasol_struct.indices{j}(i),1}.mean-mean);
         if distance(i)==0
             distance(i)=NaN;
         end
     end
     [sorted indices]=sort(distance);
     pairs{j}=parasol_struct.ids{j}(indices(1:6));
 end
 
 paired_cells=[pairs{1} pairs{2}];

end