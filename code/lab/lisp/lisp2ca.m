function ca = lisp2ca(lispstr)
% LISP2CA   Crude attempt to convert Lisp string into Cell Array
% usage: ca = lisp2ca(lispstr)
%
% Operates by converting lisp '(...)' into '{...}' and then running eval.  
% Also tries to handle lisp ':x' symbols, lisp "" strings, and lisp
% comments.  Lisp t => true, lisp nil => false.
% 
% This only works up to the level demanded by current applications; if we
% want to apply it to more complicated lisp files, additional rules will
% probably be needed.
%
% Everything is regexp based rather than PEG or CFG based.  Matlab's own
% parser should handle () well (because they are converted to {}), but ""
% cannot be handled properly.  Ultimately would be nice to replace this
% framework with a PEG approach; we could consider bringing one of the
% existing Java or C++ libraries into Matlab.
%
% 2011-05 phli
%


% Replace :x with 'x'
repstr = regexprep(lispstr, '(\s|\():([^\s]*)', '$1''$2''');

% Replace "x" with 'x'
repstr = regexprep(repstr, '"([^"]*)"', '''$1''');

% Replace parens
repstr = strrep(repstr, '#(', '{');
repstr = strrep(repstr, '(',  '{');
repstr = strrep(repstr, ')',  '}');

% Replace comments
repstr = strrep(repstr, ';', '%');

% Replace t/nil
repstr = regexprep(repstr, '(\s|{)t(\s|})', '$1true$2',  'ignorecase');
repstr = regexprep(repstr, '(\s|{)nil(\s|})', '$1false$2', 'ignorecase');

ca = eval(repstr);
