function outstructs = urgb_raster_subplot(datarun, cellspec, triggers, urgbs, varargin)
% URGB_RASTER_SUBPLOT   Wrapper for RASTER_SUBPLOT for given stimuli
% usage: outstructs = urgb_raster_subplot(datarun, cellspec, triggers, urgbs, varargin)
% 
% See also: RASTER_SUBPLOT
%
% 2012 phli
%

opts = inputParser();

% Raster opts
opts.addParamValue('start', -0.1);
opts.addParamValue('stop', mean(diff(triggers)));
opts.addParamValue('hist_bin', 0.05);
rasteroptnames = {'start', 'stop', 'hist_bin'};

% Propagating raster opts
opts.addParamValue('color',           []);
opts.addParamValue('hist_color',      []);
opts.addParamValue('hist_line_width', []);
opts.addParamValue('LineStyle',       []);
opts.addParamValue('MarkerStyle',     []);
opts.addParamValue('MarkerSize',      []);
opts.addParamValue('raster_ticks',    []);
opts.addParamValue('hist_ticks',      []);
opts.addParamValue('axopts', {});
propagating_opts = {'color', 'hist_color', 'hist_line_width', 'LineStyle', 'MarkerStyle', 'MarkerSize', 'axopts'};

opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

rasteropts = keepfields(opts, rasteroptnames{:});

stim = datarun.stimulus;
raster_size = size(urgbs);
instructs = cell(raster_size);
for i = 1:raster_size(1)
    for j = 1:raster_size(2)
        this_plot_urgbs = urgbs{i,j};
        
        % Common usage is to just give raw urgb indices here, which is
        % interpreted to mean separate PSTHs on same subplot, rather than
        % concatenated triggers.
        if ~iscell(this_plot_urgbs), this_plot_urgbs = num2cell(this_plot_urgbs); end
        
        for k = 1:length(this_plot_urgbs)
            instruct = struct();
            instruct.datarun = datarun;
            instruct.cellspec = cellspec;

            instruct.opts = rasteropts;
            for propagating_opt = propagating_opts;
                opt = propagating_opt{1};
                if isempty(opts.(opt)), continue; end
                
                if numel(opts.(opt)) == 1
                    instruct.opts.(opt) = opts.(opt);
                elseif numel(opts.(opt)) == raster_size(1)
                    instruct.opts.(opt) = opts.(opt){i};
                else
                    instruct.opts.(opt) = opts.(opt){i,j}(k);
                end
                
                if iscell(instruct.opts.(opt))
                    instruct.opts.(opt) = instruct.opts.(opt){1};
                end
            end
            
            this_concat_urgbs = this_plot_urgbs{k};
            selected_triggers = [stim.urgbi{this_concat_urgbs}];
            selected_triggers = selected_triggers(selected_triggers <= length(triggers));
            instruct.triggers = triggers(selected_triggers);
            
            instructs{i,j}(k) = instruct;
        end
    end
end

outstructs = raster_subplot(instructs, unmatched);
if nargout < 1, clear outstructs; end