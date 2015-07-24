function [s1 t1 s2 t2] = get_sep_filters(p,n,ntime,basis_kspace,basis_ktime,plotFlag)

offset = 1;

s1 = p(offset+1:offset+n);
t1 = p(offset+n+1:offset+n+ntime);

offset = offset + n + ntime;

s2 = p(offset+1:offset+n);
t2 = p(offset+n+1:offset+n+ntime);

if (~isempty(basis_kspace))
    s1 = basis_kspace*s1;
    s2 = basis_kspace*s2;
end

if (~isempty(basis_ktime))
    t1 = basis_ktime*t1;
    t2 = basis_ktime*t2;
end

if (~exist('plotFlag','var'))
    plotFlag = 0;
end

% Sign is arbitrary - by convention make the spatial profiles a positive
% bump
%fprintf('this version');
if (~determine_polarity(s1,10))
    t1 = -t1;
    s1 = -s1;
end

if (~determine_polarity(s2,10))
    t2 = -t2;
    s2 = -s2;
end


% Scale is arbitrary - by conventio make the spatial profiles unit norm
t1 = t1.* norm(s1);
s1 = s1./norm(s1);

t2 = t2.* norm(s2);
s2 = s2./norm(s2);


if (plotFlag)
   
    figure;
    subplot(2,2,1);
    imagesc(reshape(s1,sqrt(length(s1)),sqrt(length(s1)))), axis image, colorbar;
    subplot(2,2,2);
    imagesc(reshape(s2,sqrt(length(s1)),sqrt(length(s1)))), axis image, colorbar;
    subplot(2,2,3);
    plot(t1);
    subplot(2,2,4);
    plot(t2);
    
end