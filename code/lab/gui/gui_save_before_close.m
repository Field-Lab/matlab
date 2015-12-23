function gui_save_before_close(handle, datanames, varargin)
% GUI_SAVE_BEFORE_CLOSE    Checks that data has been saved before close figure
% usage: gui_save_before_close(handle, datanames, opts)
%
% opts: savenames   []   Use this to specify names that the data fields
%                        should be compared to.  Otherwise defaults to 
%                        "saved_[dataname]" for each dataname.
%
% Calling this ensures that a check is run before allowing the figure
% window to close.  The check does not verify that the data has been saved
% /per se/; that would be impossible to know.  Instead the check compares
% the current data to another GUI variable representing the state of the
% data the last time it was saved.  Therefore, the coder should include
% rewriting this save variable as part of the GUI method that actually
% saves the data.
%
% NOTE: This function relies on HAS_CALLBACK, which uses
% undocumented/unstable MatLab internal functions. This may break if new
% MatLab releases change the backend.
%
% 2010-05 phli
%

opts = inputParser;
opts.addParamValue('savenames', []);
opts.parse(varargin{:});
opts = opts.Results;
savenames = opts.savenames;

fig = getfig(handle);
if ischar(datanames)
    datanames = {datanames};
end
if ~iscell(datanames)
    error('GUI_SAVE_BEFORE_CLOSE:invalid_datanames', 'Invalid datanames argument');
end


setappdata(fig, 'save_before_close_datanames', datanames);
setappdata(fig, 'save_before_close_savenames', savenames);

if ~has_callback(fig, @save_before_closereq, 'CloseRequestFcn')
    iptaddcallback(fig, 'CloseRequestFcn', @save_before_closereq);
    api = getappdata(fig, 'api');
    api.save_before_closereq = @save_before_closereq;
    setappdata(fig, 'api', api);
end



function save_before_closereq(handle, ev) %#ok<INUSD>
fig = getfig(handle);

datanames = getappdata(fig, 'save_before_close_datanames');
savenames = getappdata(fig, 'save_before_close_savenames');

unsaved = [];
for i = 1:length(datanames);
    dataname = datanames{i};
    if ~isempty(savenames)
        savename = savenames{i};
    else
        savename = ['saved_' dataname];
    end
    
    data = getappdata(fig, dataname);
    saved = getappdata(fig, savename);
    unsaved(i) = data ~= saved;
end

if any(unsaved)
    prompt = 'Unsaved data:';
    for i = 1:length(unsaved)
        if unsaved(i)
            prompt = [prompt datanames(i)];
        end
    end
    prompt = [prompt 'Close figure anyway?'];
    
    selection = questdlg(prompt, 'Close figure?', 'Yes','No','Yes');
    if isempty(selection) || strcmp(selection, 'No')
        return
    end
end

delete(fig);