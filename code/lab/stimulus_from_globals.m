function stimulus = stimulus_from_globals(globals)

stimulus = struct();
if ~globals.runTimeMovieParamsExists(), return; end

rtmp = globals.getRunTimeMovieParams();
if rtmp.height              > 0, stimulus.field_height   = rtmp.height;             end
if rtmp.width               > 0, stimulus.field_width    = rtmp.width;              end
if rtmp.interval            > 0, stimulus.interval       = rtmp.interval;           end
if rtmp.refreshPeriod       > 0, stimulus.refresh_period = rtmp.refreshPeriod;      end
if rtmp.pixelsPerStixelY    > 0, stimulus.stixel_height  = rtmp.pixelsPerStixelY;   end
if rtmp.pixelsPerStixelX    > 0, stimulus.stixel_width   = rtmp.pixelsPerStixelX;   end

% x offset and y offset we usually blow off when we write the XML, so don't try to load them