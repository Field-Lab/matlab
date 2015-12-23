function list=is_ei_neighbor(electrode,electrode_list,array)
%greschner

if nargin<3
    array=512;
end

switch array
    case 61
        array=1;
    case 512
        array=505;
    case 519
        array=1600;
end


list=zeros(length(electrode_list),1);
list(find(electrode_list==electrode))=1;
list=ei_rev(electrode,electrode_list,array,list);




% Helper function
function list=ei_rev(electrode,electrode_list,array,list)

electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(array);  

neighbors = electrodeMap.getAdjacentsTo(electrode, 1);
t1=find(ismember(electrode_list, neighbors(2:end)));

if ~isempty(t1)
    for i=1:length(t1)
        if ~list(t1(i))
            list(t1(i))=1;
            list=ei_rev(electrode_list(t1(i)),electrode_list,array,list);
        end
    end
end