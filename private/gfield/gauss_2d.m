function f = gauss_2d(params,center,weights,locations)
%DeVries
x=locations(:,1);
y=locations(:,2);

h = center(1);%offset
k = center(2);%offset

sx = params(1);%sd
sy = params(2);%sd

A = params(3);%scale
w = params(4);%angle


fit_weights = A.*exp(-.5.*( ( ((x-h).*cos(w)-(y-k).*sin(w))./sx ).^2 + ( ((y-k).*cos(w)+(x-h).*sin(w))./sy ).^2 ) );


f = sqrt(sum((fit_weights-weights).^2));
if sx > 25 || sy > 25
    f = f + 1000;
end

if sx < 2 || sy < 2
    f = f +1000;
end

