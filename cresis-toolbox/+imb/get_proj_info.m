function [proj] = get_proj_info(loc)
% [proj] = get_proj_info(loc)

physical_constants;

switch loc
  case 'arctic'
    proj = arctic_proj;
  case 'antarctic'
    proj = antarctic_proj;
end
