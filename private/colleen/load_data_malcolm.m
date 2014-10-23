function datarun = load_data_malcolm(dataset)
    
    if(dataset==1)
        datarun=load_data('2007-03-27-1/data011-nwpca');
        datarun=load_sta(datarun);
        datarun=load_neurons(datarun);
        datarun=load_params(datarun);
        datarun=set_polarities(datarun);
    elseif(dataset==2)
        datarun=load_data('2007-08-24-4/data001-nwpca');
        datarun=load_sta(datarun);
        datarun=load_neurons(datarun);
        datarun=load_params(datarun);
        datarun=set_polarities(datarun);
    end

end

