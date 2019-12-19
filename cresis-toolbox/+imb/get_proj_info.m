function [proj] = get_proj_info(loc)
% [proj] = get_proj_info(loc)

standard_projections;

switch loc
  case 'arctic'
    proj = arctic_proj;
  case 'antarctic'
    proj = antarctic_proj;
end
