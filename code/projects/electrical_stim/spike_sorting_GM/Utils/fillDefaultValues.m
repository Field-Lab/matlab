function input = fillDefaultValues(input)


input.tracesInfo.J = size(input.tracesInfo.data,1);
input.tracesInfo.E = length(input.tracesInfo.recElecs);
input.tracesInfo.T = size(input.tracesInfo.data{1,1},2);

for j=1:input.tracesInfo.J
    
    I(j)           = size(input.tracesInfo.data{j,1},1);
end

input.tracesInfo.I = I;


input.neuronInfo.nNeurons = length(input.neuronInfo.neuronIds);

try(input.params.degPolRule)
    
catch
    input.params.degPolRule=[1 2 8 17 20 25];
end

try(input.params.Tdivide)
    
catch
    input.params.Tdivide=[0 floor(input.tracesInfo.T/2) input.tracesInfo.T];
end


    