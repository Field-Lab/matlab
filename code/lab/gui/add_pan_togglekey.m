function add_pan_togglekey(h, keys)
f = getfig(h);

iptaddcallback(h, 'WindowKeyPressFcn', @(han, ev) panon(f, ev, keys));
  
function panon(f, ev, keys)
switch ev.Key
    case keys
        pan(f, 'on');

        % Undocumented Matlab!
        %
        % Unfortunately, when pan goes on, Matlab gets really dumb and
        % disables the normal event listeners.  So we need this incantation
        % to reverse that and then add a panoff listener every single time.
        %
        % This also means the panoff listener isn't there if pan is turned
        % on one of the normal ways (toolbar button or with PAN function
        % directly).  You only get the panoff key listener if you turned
        % pan on with the key listener here.  Annoying but I can live with
        % it...
        mm = uigetmodemanager(f);
        set(mm.WindowListenerHandles,'Enable','off');
        
        set(f, 'WindowKeyPressFcn', @(han, ev) panoff(f, ev, keys));
end

function panoff(f, ev, keys)
switch ev.Key
    case keys
        pan(f, 'off');
end