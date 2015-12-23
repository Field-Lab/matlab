function data = waitfordata(f, datanames)
% WAITFORDATA   Block until figure closes, grab data from figure as it closes
% usage: data = waitfordata(f, datanames)
%
% DATANAMES is a cell array of string names for APPDATA in the figure F.
%
% Output DATA is a cell array holding the figure's APPDATA for each
% DATANAME
%
% If a naked string is passed, this is automatically boxed into a 1-element
% cell array, but note that the DATA output still comes out as a cell
% array (no unboxing).
%
% Imagine you have a GUI that saves some data you want in an APPDATA
% field named 'points'.  You can get that field using: 
%   getappdata(f, 'points')
% but this gives you no guarantee that the user is actually done with the 
% GUI; you don't know if POINTS is in a finished state.
%
% If you instead run:
%   data = waitfordata(f, {'points'});
%   points = data{1};
% you now block execution until the GUI is closed, but grab the points data
% from the GUI just as it is closing.  You have a guarantee that points
% will be set with the data from after the user finished with the GUI.

% Version 0.5
% Peter H. Li 12-Sep-2012
% As required by MatLab Central FileExchange, licensed under the FreeBSD License


% Box datanames if necessary
if ischar(datanames), datanames = {datanames}; end

% Create invisible signal checkbox
cb = uicontrol(f, 'Style',   'checkbox', ...
                  'Visible', 'off',      ...
                  'Value',   0,          ...
                  'Tag',     'waitfordata_checkbox');

% Intercept old closereq, cause it to signal the checkbox instead
oldclosereq = get(f, 'CloseRequestFcn');
set(f, 'CloseRequestFcn', @(h,e)(set(cb, 'Value', 1)));

% Wait for the checkbox to reactivate us
waitfor(cb, 'Value', 1);

% After the new closereq reactivates us; get data
data = cell(1,length(datanames));
for i = 1:length(datanames)
    data{i} = getappdata(f, datanames{i});
end

% Now run original closereqs.
if ischar(oldclosereq)
    
    % Not sure how this works, but somehow the default closereq takes no 
    % args, while any custom callback must take two.  So we must handle 
    % closereq specially.  I think handling this one special case this dumb 
    % way works, but for more industrial strength solution we could try to 
    % check the arity of the function handle...
    if strcmp(oldclosereq, 'closereq')
        closereq();
        return
    end
    
    oldclosereq = str2func(oldclosereq);
end
oldclosereq(f,[]);