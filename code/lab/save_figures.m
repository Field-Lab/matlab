function save_figures(varargin)
% SAVE_FIGURES    Save open figures
%
% usage: save_figures([opts])
%
% opts: hidden  Whether to include figures with HandleVisibility set to 
%               off, as with figures that have had LOCK called on them.  
%               Defaults to true.
%
%       format  What formats each figure should be saved as.  Passed 
%               directly to SAVEAS.  Defaults to {'fig', 'png'}.
%
% 2010-02 phli
%

opts = inputParser;
opts.addParamValue('hidden', true);
opts.addParamValue('format', {'fig', 'png'})
opts.parse(varargin{:});
opts = opts.Results;

if ischar(opts.format)
    opts.format = {opts.format};
end


orig_hidden_handle_setting = get(0, 'ShowHiddenHandles');
if opts.hidden
    set(0, 'ShowHiddenHandles', 'on');
end

try
    fhandles = get(0, 'Children');
    for i = fhandles(:)'
        filename = get_savename(i);        
        for j = 1:numel(opts.format)
            disp(['Saving figure ' num2str(i) ' as ' opts.format{j} ' to: ' filename]);
            saveas(i, filename, opts.format{j});
        end
    end
    
catch e
    do_ensure(orig_hidden_handle_setting)
    rethrow(e);
end

do_ensure(orig_hidden_handle_setting);


function do_ensure(orig_hidden_handle_setting)
set(0, 'ShowHiddenHandles', orig_hidden_handle_setting);