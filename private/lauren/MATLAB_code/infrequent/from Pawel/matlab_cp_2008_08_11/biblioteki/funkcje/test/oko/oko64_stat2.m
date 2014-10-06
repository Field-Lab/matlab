function dane=oko64_extract(name,nrchns);

name
namesize=size(name)
%name_data(1,1:namesize)=name;
%name_data(1,(namesize+1):(namesize+8))='_data.dat'


%+'_data.dat';
%name_clk=name+'_clk.dat';
%name_trig=name+'_trig.dat';

%okreslenie uzytecznego zakresu danych - wykrycie obu trigerow

