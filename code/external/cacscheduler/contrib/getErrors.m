function ers = getErrors(job, silent)
% GETERRORS(job) - Pretty print or return the errors associated with a job.
%  All task errors are collected and printed.  If silent is set to true,
%  no printing occurs and instead the bundled errors are returned. silent
%  defaults to false.
%
% EXAMPLES:
%  Standard use, pretty pring all errors from a job. 
%   getErrors(job);
%  Retrieve the errors without printing in a cell array.
%   ers = getErrors(job,true);
%
% See also getOutput, LittleJohn
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

if nargin < 2
    silent = false;
end
alltasks = get(job,'Tasks');
ers = get(alltasks,'Error');
if ~silent
    if length(ers) == 1
        if strcmp(ers.identifier,'')
            fprintf('No error!\n');
        else
            fprintf('MException! %s - %s\n', ers.identifier,ers.message);
            fprintf('Stack:\n');
            for i=1:length(ers.stack)
                fprintf('\t%s line %i\n', ers.stack(i).file,ers.stack(i).line);
            end
        end
    else
        for i = 1:length(ers)
            if strcmp(ers{i}.identifier,'')
                fprintf('Task %i - No error!\n', i);
            else
                fprintf('MException! %s - %s\n', ers{i}.identifier,ers{i}.message);
                fprintf('Stack:\n');
                for j=1:length(ers{i}.stack)
                    fprintf('\t%s line %i\n', ers{i}.stack(j).file,ers{i}.stack(j).line);
                end
            end
        end
    end
end
