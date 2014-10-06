%RawDataPath='J:\2010-09-14-0\data002min009';       
MovieFilePath='E:\data\2010-09-14-0\movie002';
patternID = 62;
electrodeID = PrimaryElectrode;

movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,patternID);
WhichMovies=[5:8];
%break
plotIndex=1;
figure(1)
for r=5%1:size(movieIDList,2);

[electrodeRawData, patternRepeatData] = NS512_GetPatternPerElectrodeData(movieIDList(r), patternID, electrodeID, full_path);

   patternInmovieRep = unique(patternRepeatData);
    for i=1:size(patternInmovieRep,1) 
        patternInmovieRepIndexes = find(patternRepeatData==patternInmovieRep(i));
        subplot(size(movieIDList,2),2,plotIndex);
        hold on
        for j=patternInmovieRepIndexes(1):patternInmovieRepIndexes(end)
            plot(electrodeRawData(j,:));
        end
        hold off
        plotIndex=plotIndex+1;
    end
end
figure(11)
plot(electrodeRawData')
break
print('-dtiff','-r300',['el' int2str(electrodeID) 'p' int2str(patternID) 'r' int2str(patternInmovieRep(i)) '.tif'])