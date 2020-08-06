function [surface, debug] = trws2_CT_perm(image, at_slope, at_weight, max_loops, ct_bounds, ft_bounds_top, ft_bounds_bottom)

% TODO[reece]: Set noise floor with echo_norm
% Set values outside bounds to noise floor
% pass into regular TRWS2 (or modify to skip columns entirely outside bounds)
% Crop output from TRWS2 to bounds

  NOISE_FLOOR = -40;
  PADEN_SOL = true;

  % Remove data entirely outside ft_bounds
  min_bound = min(ft_bounds_top(:)) + 1;
  ft_bounds_top = ft_bounds_top - min_bound + 1;
  max_bound = max(ft_bounds_bottom(:)) + 1;
  ft_bounds_bottom = ft_bounds_bottom - min_bound + 1;
  new_image = image(min_bound:max_bound, :, :);
  
  [nt, nsv, nx] = size(new_image);
  new_image = echo_norm(new_image, struct('scale', [NOISE_FLOOR 90]));

  for rline = 1:nx
    for doa_bin = 1:nsv
      new_image(1:nt < ft_bounds_top(doa_bin, rline), doa_bin, rline) = NOISE_FLOOR;
      new_image(1:nt > ft_bounds_bottom(doa_bin, rline), doa_bin, rline) = NOISE_FLOOR;
    end
  end 

  new_image = permute(new_image, [2 1 3]);
  ct_slope = zeros(size(new_image, 2), size(new_image, 3));
  ct_weight = ones(1, size(new_image, 2))*at_weight(1);

  if PADEN_SOL
    [surface, debug] = tomo.trws2_surf_bounds(single(new_image), single(at_slope), single(at_weight), single(ct_slope), single(ct_weight), uint32(max_loops), uint32(ct_bounds), uint32(ft_bounds_top), uint32(ft_bounds_bottom));
    debug = permute(debug, [2 1 3]);
    debug = [nan(min_bound - 1, nsv, nx); debug; nan(size(image, 1) - max_bound, nsv, nx)];
  else
    surface = tomo.trws2(single(new_image), single(at_slope), single(at_weight), single(ct_slope), single(ct_weight), uint32(max_loops), uint32(ct_bounds));
    surface = single(surface);
  end
  
  surface = [nan(min_bound - 1, nx); surface; nan(size(image, 1) - max_bound, nx)];
  
  for rline = 1:nx
    for rbin = 1:nt
        doa_bin = surface(rbin, rline);
        if isnan(doa_bin)
          continue;
        end
        if rbin < ft_bounds_top(doa_bin, rline) + 1 || rbin > ft_bounds_bottom(doa_bin, rline) + 1
            surface(rbin, rline) = NaN;
        end
    end
  end
end