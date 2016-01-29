if 1
    if 1
        clear datarun
 
        if 1
            
           % exp_nm = '2012-09-27-3';  
         %  exp_nm = '2012-08-09-3'; dn_mas = sprintf('data002'); dn_slvBW = sprintf('data006'); dn_slvNSEM = sprintf('data005');
           exp_nm = '2012-09-21-1'; dn_mas = sprintf('data004'); dn_slvBW = sprintf('data006'); dn_slvNSEM = sprintf('data005');
           exp_nm = '2012-09-24-2'; dn_mas = sprintf('data002'); dn_slvBW = sprintf('data005'); dn_slvNSEM = sprintf('data001');
           exp_nm = '2012-09-13-1'; dn_mas = sprintf('data004'); dn_slvBW = sprintf('data007'); dn_slvNSEM = sprintf('data003');
           dn_slv = dn_slvNSEM;  dn_slv = dn_slvBW;
            com_dir = sprintf('/Akheitman/Analysis/%s',exp_nm); % COMMON DIRECTORY FOR MASTER SLAVE
        %    dn_mas = sprintf('data003');  % CHOOSE APPROPRIATE BOOKEND
           % dn_slv = sprintf('data002');

            datarun{1}.piece.array_id=1501;           
            datarun{1}.names.rrs_params_path=sprintf('%s/%s/%s.params',com_dir,dn_mas,dn_mas);
            datarun{1}.names.rrs_neurons_path=sprintf('%s/%s/%s.neurons',com_dir,dn_mas,dn_mas);
            datarun{1}.names.rrs_ei_path=sprintf('%s/%s/%s.ei',com_dir,dn_mas,dn_mas);
            datarun{1}.names.shortname  =dn_mas;
            %datarun{1}.piece.array_id=1501;
           % datarun{1}.names.rrs_params_path='/snle/lab/Experiments/Array/Analysis/2011-08-04-7/data000/data000.params';
           % datarun{1}.names.rrs_neurons_path='/snle/lab/Experiments/Array/Analysis/2011-08-04-7/data000/data000.neurons';
           % datarun{1}.names.rrs_ei_path='/snle/lab/Experiments/Array/Analysis/2011-08-04-7/data000/data000.ei';
           % datarun{1}.piece.array_id=504;
           
           
            datarun{2}.names.rrs_neurons_path=sprintf('%s/%s/%s.neurons',com_dir,dn_slv,dn_slv);
            datarun{2}.names.rrs_ei_path=sprintf('%s/%s/%s.ei',com_dir,dn_slv,dn_slv);
            datarun{2}.piece.array_id=1501;
            datarun{2}.names.shortname  =dn_slv;

           

            %datarun{2}.names.rrs_neurons_path='/snle/lab/Experiments/Array/Analysis/2011-07-14-4/data006-007-nwpca/data006-007-nwpca.neurons';
            %datarun{2}.names.rrs_ei_path='/snle/lab/Experiments/Array/Analysis/2011-07-14-4/data006-007-nwpca/data006-007-nwpca.ei';
            %datarun{2}.piece.array_id=1501;

            %save_path='/snle/lab/Experiments/Array/Analysis/2012-07-26-0/data002-nwpca/data002-nwpca/';

            save_path=sprintf('%s/%s/',com_dir,dn_slv);
        end    

        opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_ei',1);
        datarun=load_data(datarun,opt);
    end
    
    map_ei_classification_txtAH(datarun{1}, datarun{2},'master_cell_type',{1 2 3 4 5},'classification_path', save_path,...
        'corr_threshold',.7);
end


