function Wcp = wc_preprocess_javamovie_rawmovie(Wc, stimulus)
% Pre-permute Wc to avoid permuting every frame...
%
% Format apparently the same for Java movie and stimulus raw movies (not so
% surprising).
%
% 2012-09 phli
%

ind = 1:size(Wc,1);
ind = reshape(ind, stimulus.field_height, stimulus.field_width, 3);
ind = permute(ind, [3 2 1]);
ind = reshape(ind, [], 1);
Wcp = Wc(ind,:);