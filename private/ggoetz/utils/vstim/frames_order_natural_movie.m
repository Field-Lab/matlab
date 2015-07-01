function fo = frames_order_natural_movie(frames_a, frames_b, repeats)
%FRAMES_ORDER_NATURAL_MOVIE Computes the order in which the frames of a
%natural movie were displayed.
%
%  FO = FRAMES_ORDER_NATURAL_MOVIE(FA, FB, R) calculates how the frames of
%  a natural movie were interleaved before they were displayed. The movie
%  shown should have consisted of two chunks A and B of length FA and FB
%  frames repeatedly. Chunk A should have been shown repeatedly (to create
%  rasters) while chunk B should have been varied. The sequence [chunkA
%  chunkB] was repeated R times.
%  Then, the frames from the original movie file were shown in the order
%  specified by the returned value FO.
%
%  Indices are 1-based, so that if a 0 is returned it means no frame (i.e. 
%  a grey screen) was shown.

start_frames_b = frames_a + 1;

fo = [];
for i = 1:repeats
    fo = [fo 1:frames_a zeros(1, 100) (1:frames_b)+start_frames_b zeros(1, 100)]; %#ok<AGROW>
    start_frames_b = start_frames_b + frames_b;
end

end % frames_order_natural_movie