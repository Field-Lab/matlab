function NumberOfMovies = NS512_GetNumberOfMovies_GlobalPath(full_path)

%check source directory
if (full_path(length(full_path))=='\') separator=''; else separator='\'; end
files = dir([full_path separator 'data*.bin']);
%if(size(files, 1)==0) error(['There are no files to process in ' full_path]); end

if (full_path(length(full_path))=='\') shiftN=1; else shiftN=0; end
inputNumber = full_path(length(full_path)-2-shiftN:length(full_path)-shiftN)

cd (full_path)
cd ..

%fid0=fopen(['movie' inputNumber],'r','b')
fid0=fopen(['movie' inputNumber],'r','b');
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies;
fclose(fid0);

end

