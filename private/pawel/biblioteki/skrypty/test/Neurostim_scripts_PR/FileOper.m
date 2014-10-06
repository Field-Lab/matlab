%Universal, custom file operation script
CurrentPath =  'D:\Home\Rydygier\NEURO\analysis\retina\2009-11-27-0\data001\April2010\figures\50us_diff_met5met6';
cd(CurrentPath);
ListOfFilesStructure = dir;
FirstFileNumber = 3;
FileName = [];
p = [];
m = [];
el = [];

for i=FirstFileNumber:length(ListOfFilesStructure)
    FileName = ListOfFilesStructure(i).name;
    
    pIndex = strfind(FileName, 'p');
    mIndex = strfind(FileName, 'm');
    eIndex = strfind(FileName, 'e');
    dotIndex = strfind(FileName, '.');
    
        p = [p str2num(FileName(pIndex+1:mIndex-2))];
        m = [m str2num(FileName(mIndex+1:eIndex-2))];
        el = [el str2num(FileName(eIndex+2:dotIndex-1))];
end
p = p';
m = m';
el = el';