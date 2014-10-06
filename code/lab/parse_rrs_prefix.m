function parsed = parse_rrs_prefix(rp)
% PARSE_RRS_PREFIX      Get piece/run/etc. from rrs_prefix
% usage: parsed = parse_rrs_prefix(rp)
%
% Gives back a struct with elements year, month, day, piece, run, etc.
%
% 2011-07 phli
%

parsed = struct([]);
rgxp = ['\' filesep '(\d{4})-(\d{2})-(\d{2})-(\d).*\' filesep '(data|d)(\d{2,3})(.*)\' filesep];
tokens = regexp(rp, rgxp, 'tokens');
if length(tokens) ~= 1, return; end

tokens = tokens{1};
parsed = struct();
parsed.year  = tokens{1};
parsed.month = tokens{2};
parsed.day   = tokens{3};
parsed.piece = tokens{4};
parsed.piece_fullname = [parsed.year '-' parsed.month '-' parsed.day '-' parsed.piece];
parsed.piece_num  = str2double(parsed.piece);

parsed.run      = tokens{6};
parsed.run_num  = str2double(parsed.run);
parsed.run_name = [tokens{5} tokens{6}];
parsed.run_full_name = [parsed.run_name tokens{7}];