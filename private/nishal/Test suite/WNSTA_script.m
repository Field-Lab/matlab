

WN_datafile = '2015-03-09-3/streamed/data008/data008';
WN_datafile_short='2015-03-09-3/streamed/data008/data008';
movie_xml = 'RGB-8-4-0.48-11111';
stim_len=900;% in seconds

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


for cellTypeId=1 % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 




for ref_cell_number=1
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    
 
  
  xx=datarun.cell_types{cellTypeId}.name;
  xx(xx==' ') = '';
  dsave=sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data008/CellType_%s/rk2',xx);
  Mask=zeros(80,40);

 [WNSTA,STC,Mask]=STA_STC_from_WNrun({cellID}, WN_datafile, movie_xml, 900,dsave,Mask)

    InterestingCell_vis_id(ref_cell_number)

end


end

%%
WN_datafile = '2014-08-20-0/data000_wrong_nps/data000';
movie_xml = 'RGB-4-1-0.48-11111';
stim_len=1800;% in seconds
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
cellID=2673;
dsave=sprintf('/Volumes/Lab/Users/bhaishahster/STC_computation');
Mask=zeros(80,40);
[WNSTA,STC,Mask]=STA_STC_from_WNrun({cellID}, WN_datafile, movie_xml, stim_len,dsave,Mask)
InterestingCell_vis_id(ref_cell_number)

%% 
WN_datafile = '2015-03-09-3/streamed/data008/data008';
movie_xml = 'RGB-8-4-0.48-11111';
stim_len=900;% in seconds
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
cellID=3121;
dsave=sprintf('/Volumes/Lab/Users/bhaishahster/STC_computation');
Mask=zeros(80,40);
Mask(19-10:19+10,27-10:27+10)=1;
[WNSTA,dimensions,Mask]=STA_STC_from_WNrun({cellID}, WN_datafile, movie_xml, stim_len,dsave,Mask)
InterestingCell_vis_id(ref_cell_number)

