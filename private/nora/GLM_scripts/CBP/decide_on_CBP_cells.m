for i=2:4
    [piece,cells,expname]=cell_list(i,'shortlist');
    number_of_cells=length(cells);
    [exp_info]=experimentinfoNSEM(piece);
    datarun=load_data([piece '/' exp_info.dr.mas], struct('load_neurons',true));
    rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(['/Volumes/Data/' piece '/' exp_info.dr.mas]);
    figure(i);
    for k=1:number_of_cells
        subplot(4,3,k)
        channel=datarun.channels(datarun.cell_ids==cells{k});
        data=rdf.getData(channel,1,20000);
        plot(data)
        title(num2str(cells{k}))
        ylim([-200 400])
    end
end

CBPcells{1}.piece='2012-08-09-3';
CBPcells{1}.ON=1772;
CBPcells{1}.OFF=1471;
CBPcells{1}.WN='data006';
CBPcells{1}.NSEM='data005';

% CBPcells{2}.piece='2012-09-27-3';
% CBPcells{2}.ON=1909;
% CBPcells{2}.OFF=1;
% CBPcells{2}.WN='data005';
% CBPcells{2}.NSEM='data002';

CBPcells{3}.piece='2013-08-19-6';
CBPcells{3}.OFF=1328;
CBPcells{3}.ON=3167;
CBPcells{3}.WN='data003';
CBPcells{3}.NSEM='data001';
CBPcells{3}.ONcoupling=[3396,3350,3137,2855,3694,3306];

CBPcells{2}.piece='2013-10-10-0';
CBPcells{2}.OFF=346;
CBPcells{2}.ON=32;
CBPcells{2}.WN='data005';
CBPcells{2}.NSEM='data001';

save('CBPcells.mat','CBPcells')