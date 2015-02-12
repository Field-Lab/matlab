function strct = keepfields(strct, varargin)
% keepfields    Remove all but the specified fields from the struct

force = false;
if strcmp(varargin{end}, '-force')
    force = true;
    varargin = varargin(1:end-1);
end

if iscell(varargin{1})
    keptfields = varargin{1};
else
    keptfields = varargin;
end


if ~force
    notfound = setdiff(keptfields, fieldnames(strct));
    if ~isempty(notfound)
        error('Could not find all the specified fields to keep; was there a typo?');
    end
end

rmfields = setdiff(fieldnames(strct), keptfields);
strct = rmfield(strct, rmfields);
