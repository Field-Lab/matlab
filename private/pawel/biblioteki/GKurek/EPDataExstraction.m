clear
clc

full_path='J:\2010-09-14-0\data002min009';       
patternID = 193;
electrodeID = 12;

%movieIDList = NS512_GetMoviesWithPatternID(full_path,patternID);
SelectedMovies = NS512_GetMoviesWithPatternIDNew(full_path,patternID);
movieIDList = [1 9 17 25 33];

plotIndex=1;
for r=1:size(movieIDList,2)

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

print('-dtiff','-r300',['el' int2str(electrodeID) 'p' int2str(patternID) 'r' int2str(patternInmovieRep(i)) '.tif'])