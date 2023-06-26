function  eps_r = perm_profile(h,time,mode,er)

if strcmpi(mode,'constant')
  air_ice_idx  = find(time >= h,1);
  if isempty(air_ice_idx)
    % all air
    eps_r        = ones(length(time),1)*1;
  else
    eps_r        = [ones(air_ice_idx,1)*1; ones(length(time)-air_ice_idx,1)*er];
  end
elseif strcmpi(mode,'real')
  % FILL THIS IN
elseif strcmpi(mode,'complex')
  % FILL THIS IN
else
  error('Unsupported:Mode','Unsupported mode')
end

return