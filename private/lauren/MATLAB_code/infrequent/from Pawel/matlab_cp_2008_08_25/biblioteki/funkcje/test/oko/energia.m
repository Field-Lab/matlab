function y=energia(signal,offset);
%Zwraca usredniona energie sygnalu przyjmujac parametr "offset" jako DC, czyli: 
%- odejumje offset od calego sygnalu "signal";
%- sumuje probki z kwadratem;
%- dzieli przez ilosc probek.
%sgsdgsg
s=signal-offset;
y=sum(s.^2)/length(s);