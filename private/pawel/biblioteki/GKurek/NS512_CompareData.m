clear
%clc

%full_path='E:\ne\data003';          
%full_pathTTX='E:\ne\data005';        
%output_path = 'E:\ne\data000';     
%vision_path = 'E:\ne\Vision.jar'; 

%javaaddpath(vision_path);  

full_path='G:\2010-09-14-0\data002';          
full_pathTTX='G:\2010-09-14-0\data009\';        
output_path = 'D:\Home\Pawel\analysis\2010-09-14-0\data002min009';      

%full_path='E:\ne\data003';          
%full_pathTTX='E:\ne\data005';        
%output_path = 'E:\ne\data000';   

StartPosition=69100000;      %<---
AmountOfData=400000;         %<---
sampleChannel = 40;

RawDataIN= NS512_GetData(full_path, sampleChannel, StartPosition, AmountOfData);
RawDataOUT= NS512_GetData(output_path, sampleChannel, StartPosition, AmountOfData);
%RawDataTTX= NS512_GetData(full_pathTTX, sampleChannel, StartPosition, AmountOfData);
    
figure(10)
clf
g=plot(RawDataIN);
set(g(1),'Displayname','RawData');
hold on
g=plot(RawDataOUT,'r');
set(g(1),'Displayname','RawData-TTX');

%g=plot(RawDataTTX,'g');
%set(g(1),'Displayname','TTX');

legend('Location','NorthEast')