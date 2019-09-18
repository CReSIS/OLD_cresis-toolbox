function [varargout] = nanmin(varargin)
% [varargout] = nanmin(varargin)

[varargout{1:nargout}] = min(varargin{:});
