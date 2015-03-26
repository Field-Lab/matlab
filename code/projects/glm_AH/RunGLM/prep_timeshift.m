% AKHeitman 2014-05-13
% rotate the X_bin by shifts
% need to comment more fully !!

% X_bin = 
%[a b c d e f g  
% h i j k l m n]
% shifts = 
% [2,3]

% ANswer  X_bin_shift = 
% [ f g a b c d e 
%   m n h i j k l
%   e f g a b c d 
%   l m n h i j k ]



function X_bin_shift = prep_timeshift(X_bin,shifts)

bins = size(X_bin,2);
dim  = size(X_bin,1);

X_bin_shift = zeros(dim*length(shifts) , bins);

for i_shift = 1:length(shifts)
    step = shifts(i_shift);
    index = (i_shift-1)*dim +  [1:dim] ;
    X_bin_shift(index,:) = circshift(X_bin , [0 step]);
end

end
    