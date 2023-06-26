function [status,result] = robust_save(fn,varargin)

status = -1;
while status ~= 0
  try
    evalin('caller',[sprintf('save(''%s''',fn),sprintf(',''%s''',varargin{:}),');']);
  catch ME
    ME.getReport
    fn
    warning('save call failed');
    keyboard
  end
  try
    load(fn)
    status = 0;
  catch ME
    ME.getReport
    fn
    warning('load call failed');
    keyboard
  end
end
