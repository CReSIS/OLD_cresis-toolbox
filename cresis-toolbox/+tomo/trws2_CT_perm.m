function [surface] = trws2_CT_perm(image, at_slope, at_weight, max_loops, ct_bounds, ft_bounds_top, ft_bounds_bottom)

% TODO[reece]: Set noise floor with echo_norm
% Set values outside bounds to noise floor
% pass into regular TRWS2 (or modify to skip columns entirely outside bounds)
% Crop output from TRWS2 to bounds

  [nt, nsv, nx] = size(image);

  NOISE_FLOOR = -40;

  image = echo_norm(image, struct('scale', [NOISE_FLOOR 90]));

  for rline = 1:nx
    for doa_bin = 1:nsv
      image(1:nt < ft_bounds_top(doa_bin, rline), doa_bin, rline) = NOISE_FLOOR;
      image(1:nt > ft_bounds_bottom(doa_bin, rline), doa_bin, rline) = NOISE_FLOOR;
    end
  end 

  new_image = permute(image, [2 1 3]);
  ct_slope = zeros(size(new_image, 2), size(new_image, 3));
  ct_weight = ones(1, size(new_image, 2));
  surface = tomo.trws2_surf_bounds(single(new_image), single(at_slope), single(at_weight), single(ct_slope), single(ct_weight), uint32(max_loops), uint32(ct_bounds), uint32(ft_bounds_top), uint32(ft_bounds_bottom));
   
%   for rline = 1:nx
%     for rbin = 1:nt
%         doa_bin = surface(rbin, rline);
%         if rbin < ft_bounds_top(doa_bin, rline) || rbin > ft_bounds_bottom(doa_bin, rline)
%             surface(rbin, rline) = NaN;
%         end
%     end
%   end
end