function [D]=dac(amplitude_range, current_range, relative_amplitude, relative_amplitude_max)
% amplitude_range - amplituda sygna³u max
% current_range - zakres
% relative_amplitude - wzglêdna amplituda
% relative_amplitude_max - wzglêdna amplituda maksymalna
a = 127*(amplitude_range/current_range);
b=a*abs(relative_amplitude/relative_amplitude_max);
c=b/relative_amplitude;
d=abs(int8(c));
e=double(abs(d*relative_amplitude));
f=b-e;
if f<=(abs(relative_amplitude/relative_amplitude_max)/2);
   d=d*abs(relative_amplitude);
else
    d=(d+1)*abs(relative_amplitude);
end
g=de2bi(d,7);
h=fliplr(g);
if relative_amplitude <0
    A=[0,h];
else
    A=[1,h];
end
D=[A(1:4)', A(5:8)'];
end




