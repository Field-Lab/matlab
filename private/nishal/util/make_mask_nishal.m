

function make_mask_nishal(datafile,cell_params,save_file)


%% Load cells and STAs
global_vars
% global_vars2
%datafile = '2013-10-10-0/data000';
type_name= cell(1,1);
type_name{1}=cell_params.type_name_inp;

if(~strcmp(datafile,'load_from_cell_params'))
    datarun=load_data(datafile)
    datarun=load_sta(datarun)
    datarun=load_params(datarun)
    
    %get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
    %STAs!
    if(strcmp(type_name{1},'userCellList'))
        idx=[1:length(datarun.cell_ids)];
        idx_list=[];
        for icell_check=1:length(cell_params.cell_list)
            idx_list=[idx_list;idx(datarun.cell_ids==cell_params.cell_list(icell_check))];
        end
        
        matlab_cell_ids=idx_list;
        clear idx_list
    else
        matlab_cell_ids=get_cell_indices(datarun,type_name);
    end
    stas=datarun.stas.stas(matlab_cell_ids);
else
    stas=cell_params.stas;
    matlab_cell_ids=1;
end

% Load STAs

stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:size(stas{1},4)
        %st_temp(:,:,:,itime)=mean(repmat(stas{icell}(:,:,3,end-itime+1),[1,1,3,1]),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)';
    end
    %sprintf('Only blue gun selected!!!!!')
    stas_new{icell}=st_temp;
end
stas=stas_new;
stas_orig=stas_new;


[stas_clipped,totalMaskAccept,CellMasks]= clipSTAs(stas,cell_params);
mov_params.totalMaskAccept=totalMaskAccept;
cell_params.CellMasks=CellMasks;

if(isfield(cell_params,'use_fits')==1)
    if(cell_params.use_fits==1)
        addpath(genpath('~/Nishal/matlab/code'));
        addpath(genpath('~/Nishal/matlab/private/nishal/fwdfittingfunctions'));
        
        fit_info=cell(length(stas),1);
        display('Starting STA fitting');
        parfor issta=1:length(stas)
            issta
            fit_info{issta} = fit_sta(stas_new{issta});
        end
        
        
        full_fit=cell(1,1);
        cellSelected=zeros(length(stas),1);
        icnt=0;
        for issta=1:length(stas)
            if(fit_info{issta}~=[])
                icnt=icnt+1;
                full_fit{icnt} = sta_fit_function(fit_info{issta}.initial_params);
                cellSelected(issta)=1;
            else
                cellSelected(issta)=0;
                display(sprintf('Cell Removed %d',issta));
            end
            
        end
        stas=full_fit;
        mov_params.totalMaskAccept=ones(size(stas{1},1),size(stas{1},2));
        display('Using STA fits')
        
    
            [stas_clipped,totalMaskAccept,CellMasks]= clipSTAs(stas,cell_params);
            mov_params.totalMaskAccept=totalMaskAccept;
            cell_params.CellMasks=CellMasks;
    
    end
    
    if(cell_params.use_fits==2)
        
            [stas,totalMaskAccept,CellMasks]= clipSTAs(stas,cell_params);
            mov_params.totalMaskAccept=totalMaskAccept;
            cell_params.CellMasks=CellMasks;
            stas_clipped=stas;
      
    end
end
    
cells = datarun.cell_ids(matlab_cell_ids);
mask = zeros(size(totalMaskAccept,1),size(totalMaskAccept,2),length(cells));
for icell=1:length(cells)
    mask(:,:,icell) = CellMasks{icell};
end

final_mask = double(sum(mask,3)>0);

cells =num2cell(cells);
 save(save_file,'final_mask', 'mask','cells');
end

    