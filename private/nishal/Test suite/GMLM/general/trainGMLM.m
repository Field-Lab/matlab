function GLMparams=trainGMLM(response,movie,numFilters)

% response at higher resolution and movie at lower resolution

binPerFrame=10;
HistLen=120; % in terms of bins
Histbasis=10;

Filtdim1=size(movie,1);
Filtdim2=size(movie,2);
Filtlen=30;
stas=cell(numFilters,1);

% Make movie form RGB to BW
movie=squeeze(mean(movie,3));

% Initialize variables
for icell=1:numFilters
stas{icell}=rand(Filtdim1,Filtdim2,1,Filtlen);
end
H=zeros(Histbasis,1);

[gradstas,gradH] = gradientGMLM(stas,movie,H,response,HistLen);

end

function [gradstas,gradH]=gradientGMLM(stas,movie,H,response,HistLen)

% find k1'x, k2'x, h.y ? 
n_cell=length(stas);
movie_time=size(movie,3);
kx=Ax(stas,mov,movie_time,n_cell);
hy = PSF(H,response,HistLen);
 
end

function hy=PSF(H,response,HistLen)
Histbasis=length(H);
RCBasis=zeros(Histbasis,HistLen);
for ibasis=1:Histbasis
RCBasis(ibasis,:)=
end
end