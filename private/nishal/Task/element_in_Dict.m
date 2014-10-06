function I=element_in_Dict(ji,j_cuei,Dict)
store = Dict(ji,:);
I=0;
for i=1:length(store)
I=I+double(store(i)==j_cuei);

end

end