% 2014-02-11 
% AkHeitman
% Hack to get around stupid malfunctioning stim_design c-code .. don't know
% why it broke!

% Used by linearfilt_grad3AH in gln_AH_20 on not by linearfilt_grad2AH 
% not used in glm_AH_19

%{
x = [1:5 ,5:-1:1];
shifts = 5;
L = stimdesign_hackAH(x,shifts);
imagesc(L); colorbar
%}

function L = stimdesign_hackAH(vector , shifts);

L = zeros( shifts, length(vector) );
for i_shift = 1 : shifts
    L(i_shift,:) = circshift ( vector , [0 , (i_shift-1) ] );
end
    
    

