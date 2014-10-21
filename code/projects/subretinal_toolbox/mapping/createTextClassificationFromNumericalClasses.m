function createTextClassificationFromNumericalClasses(neuronIDList, classesVector)
% This function outputs to the command line a Vision-formatted
% classification text which can then be copy and pasted into a vision
% classification file. 
% The parameters neuronIDList and classesVector should be two vectors
% that have exactly the same number of elements. 
%
% Parameters:
%   - neuronIDList: a n x 1 or 1 x n vector of neuron IDs.
%   - classesVector: a n x 1 or 1 x n vector of classes, represented by 
%   distinct integers. 
%
% Returns:
%   []
%
% Example:
%   If neurons 1, 2 and 3 belong to class 0 while neuron 4 belongs to class
%   1, the function should be called as follows:
%         createTextClassificationFromNumericalClasses([1 2 3 4], [0 0 0 1])

nNeurons = length(neuronIDList);
if length(classesVector)~=nNeurons
    err = MException('MATLAB:InvArgIn',...
        'Invalid input size');
    throw(err);
end

for kk=1:nNeurons
    display(sprintf('%d  All/Class%d',neuronIDList(kk), classesVector(kk)));
end