function [varargout] = nanmax(varargin)
% [varargout] = nanmax(varargin)

[varargout{1:nargout}] = max(varargin{:});
