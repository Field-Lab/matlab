% 1)
clear
clc

full_path='G:\data\2013-12-12-0-PH\data001';          
full_pathTTX='G:\data\2013-12-12-0-PH\data009';        
output_path = 'I:\analysis\slices\2013-12-12-0-PH\data001';
MovieFilePath='G:\data\2013-12-12-0-PH\movie001';

log_path = 'I:\analysis\slices\2013-12-12-0-PH\log001.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);

% 2)
clear
clc

full_path='G:\data\2013-12-12-0-PH\data003';          
full_pathTTX='G:\data\2013-12-12-0-PH\data007';        
output_path = 'I:\analysis\slices\2013-12-12-0-PH\data003';
MovieFilePath='G:\data\2013-12-12-0-PH\movie003';

log_path = 'I:\analysis\slices\2013-12-12-0-PH\log003.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);


% 3)
clear
clc

full_path='G:\data\2013-12-12-0-PH\data004';          
full_pathTTX='G:\data\2013-12-12-0-PH\data009';
output_path = 'I:\analysis\slices\2013-12-12-0-PH\data004';
MovieFilePath='G:\data\2013-12-12-0-PH\movie004';

log_path = 'I:\analysis\slices\2013-12-12-0-PH\log004.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);


% 4)
clear
clc

full_path='G:\data\2013-12-12-0-PH\data005';          
full_pathTTX='G:\data\2013-12-12-0-PH\data007';
output_path = 'I:\analysis\slices\2013-12-12-0-PH\data005';
MovieFilePath='G:\data\2013-12-12-0-PH\movie005';

log_path = 'I:\analysis\slices\2013-12-12-0-PH\log005.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);

% ------------

% 1)
clear
clc

full_path='G:\data\2013-12-12-3-PH\data001';          
full_pathTTX='G:\data\2013-12-12-3-PH\data010';        
output_path = 'I:\analysis\slices\2013-12-12-3-PH\data001';
MovieFilePath='G:\data\2013-12-12-3-PH\movie001';

log_path = 'I:\analysis\slices\2013-12-12-3-PH\log001.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);


% 2)
clear
clc

full_path='G:\data\2013-12-12-3-PH\data002';          
full_pathTTX='G:\data\2013-12-12-3-PH\data008';        
output_path = 'I:\analysis\slices\2013-12-12-3-PH\data002';
MovieFilePath='G:\data\2013-12-12-3-PH\movie002';

log_path = 'I:\analysis\slices\2013-12-12-3-PH\log002.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);

% 3)
clear
clc

full_path='G:\data\2013-12-12-3-PH\data003';          
full_pathTTX='G:\data\2013-12-12-3-PH\data009';        
output_path = 'I:\analysis\slices\2013-12-12-3-PH\data003';
MovieFilePath='G:\data\2013-12-12-3-PH\movie003';

log_path = 'I:\analysis\slices\2013-12-12-3-PH\log003.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);


% 4)
clear
clc

full_path='G:\data\2013-12-12-3-PH\data004';          
full_pathTTX='G:\data\2013-12-12-3-PH\data006';        
output_path = 'I:\analysis\slices\2013-12-12-3-PH\data004';
MovieFilePath='G:\data\2013-12-12-3-PH\movie004';

log_path = 'I:\analysis\slices\2013-12-12-3-PH\log004.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);


% 5)
clear
clc

full_path='G:\data\2013-12-12-3-PH\data005';          
full_pathTTX='G:\data\2013-12-12-3-PH\data007'; 
output_path = 'I:\analysis\slices\2013-12-12-3-PH\data005';
MovieFilePath='G:\data\2013-12-12-3-PH\movie005';

log_path = 'I:\analysis\slices\2013-12-12-3-PH\log005.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);



% -----------------

% 1)
clear
clc

full_path='G:\data\2013-12-15-0\data004';          
full_pathTTX='G:\data\2013-12-15-0\data011'; 
output_path = 'I:\analysis\slices\2013-12-15-0\data004';
MovieFilePath='G:\data\2013-12-15-0\movie004';

log_path = 'I:\analysis\slices\2013-12-15-0\log004.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);

% 2)
clear
clc

full_path='G:\data\2013-12-15-0\data005';          
full_pathTTX='G:\data\2013-12-15-0\data012'; 
output_path = 'I:\analysis\slices\2013-12-15-0\data005';
MovieFilePath='G:\data\2013-12-15-0\movie005';

log_path = 'I:\analysis\slices\2013-12-15-0\log005.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);


% 3)
clear
clc

full_path='G:\data\2013-12-15-0\data008';          
full_pathTTX='G:\data\2013-12-15-0\data011'; 
output_path = 'I:\analysis\slices\2013-12-15-0\data008';
MovieFilePath='G:\data\2013-12-15-0\movie008';

log_path = 'I:\analysis\slices\2013-12-15-0\log008.txt';

NumberOfMovies=NS512_GetNumberOfMovies(full_path)
logFid = fopen(log_path, 'w');

tic
  for MovieID=1:NumberOfMovies
    OutputData=NS512_ProcessRawDataAndTTX_v2(full_path, full_pathTTX, output_path, MovieFilePath, MovieID, NumberOfMovies, logFid);
  end
toc
fclose(logFid);