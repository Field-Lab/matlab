function M = remap_event_indices(M, i_pre, i_post)
%REMAP_EVENT_INDICES Remaps the indices associated with some events.
%
%  REMAP_EVENT_INDICES(M, I_PRE, I_POST) remaps a matrix of indices M 
%  belonging to a set of indices I_PRE to a set of indices I_POST. 
%  The indices are remapped so that if M(j,k) == I_PRE(l) then after
%  remapping M(j,k) == I_POST(l).
%
%  zero indices are kept zero, as they always mean no frame was shown.

for kk = 1:numel(M)
    if M(kk)
        M(kk) = i_post(i_pre == M(kk));
    end
end

end % remap_event_indices