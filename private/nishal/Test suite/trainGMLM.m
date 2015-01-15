function GLMparams=trainGMLM(response,movie,numFilters)

binPerFrame=10;
HistLen=12; % in terms of bins
H=zeros(HistLen,1);

Filtdim1=size(movie,1);
Filtdim2=size(movie,2);
Filtlen=30;
stas=cell(numFilters,1);

end
