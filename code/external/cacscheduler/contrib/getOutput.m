function out = getOutput(job,silent)
% GETOUTPUT(job) - Pretty print or return the command window output for all
%  tasks in the provided job object.  If silent is set to true, output is
%  returned in a cell array without any printing.  silent defaults to
%  false.
%
% EXAMPLES:
%  Standard use, pretty pring all output from a job.
%   getOutput(job);
%  Retrieve the output without printing in a cell array.
%   ers = getOutput(job,true);
%
% See also getErrors, LittleJohn
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

if nargin < 2
    silent = false;
end
alltasks = get(job,'Tasks');
out = get(alltasks,'CommandWindowOutput');

if ~silent
    if length(alltasks) == 1
        %Verify that the task actually had CaptureCommandWindowOutput On
        stat = get(alltasks(1),'CaptureCommandWindowOutput');
        if stat == 0
            fprintf('CaptureCommandWindowOutput is false! Set to true prior to job submission!\b');
        end
        fprintf('Task %d out: %s\n',1,out);
    else
        for i = 1:length(alltasks)
            stat = get(alltasks(i),'CaptureCommandWindowOutput');
            if stat == 0
                fprintf('CaptureCommandWindowOutput is false for task %i! Set to true prior to job submission!\b',i);
            end
            fprintf('Task %d out:%s\n', i,out{i});
        end
    end
end
