function [output] = prefilter(maskedMovdd_sliced,filters)
output =zeros(length(filters),size(maskedMovdd_sliced,2));

for ifilter = 1:length(filters)
output(ifilter,:) = filters{ifilter}'*maskedMovdd_sliced;    
end

end