function [status,result] = robust_system(cmd)

status = -1;
attempts = 0;
while status ~= 0
  try
    [status,result] = system(cmd);
  catch ME
    ME.getReport
    cmd
    warning('system call failed');
    keyboard
  end
  if status ~= 0
    warning('%s\n  FAILED %d %s', cmd, status, result);
    attempts = attempts + 1;
    delay_period = 3*2^(attempts-1);
    if delay_period > 120
      delay_period = 120;
    end
    fprintf('  Delaying %d seconds (%s)\n', delay_period, datestr(now));
    pause(delay_period);
    if attempts > 10
      warning('There is potentially something wrong with the system command if it fails this many times. Look into it before running "dbcont". system command is:\n  %s', cmd);
      keyboard;
    end
  end
end
