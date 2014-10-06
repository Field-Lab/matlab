function varargout = plot_xy(xy, varargin)
varargout{:} = plot(xy(:,1), xy(:,2), varargin{:});