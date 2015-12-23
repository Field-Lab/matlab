% akheitman 2014 - 03-05
% Just book keeipng .. so we sort by cell ids and then find type easily
function [ celltype , cell_sname, cell_lname]  = findcelltype(cid, datarun_cell_types)


if  ~isempty(  find(datarun_cell_types{1}.cell_ids == cid) )
    disp('here')
        celltype = 'ON-Parasol';
        celltype_sname = 'ONPar';
end
if  ~isempty(  find(datarun_cell_types{2}.cell_ids == cid) )
        celltype = 'OFF-Parasol';
        celltype_sname = 'OFFPar';
end
if  ~isempty(  find(datarun_cell_types{3}.cell_ids == cid) )
        celltype = 'ON-Midget';
        celltype_sname = 'ONMid';
end
if ~isempty(  find(datarun_cell_types{4}.cell_ids == cid) )
        celltype = 'OFF-Midget';
        celltype_sname = 'OFFMid';
end
if ~isempty(  find(datarun_cell_types{5}.cell_ids == cid) )
        celltype = 'SBC';
        celltype_sname = 'SBC';
end


cell_sname = sprintf('%s_%d', celltype_sname, cid) ;
cell_lname = sprintf('%s: Cell-%d', celltype, cid);
end

    
