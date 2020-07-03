function [surface] = trws2_CT_perm(image, at_slope, at_weight, max_loops, bounds)

  new_image = permute(image, [2 1 3]);
  ct_slope = zeros(size(new_image, 2), size(new_image, 3));
  ct_weight = ones(1, size(new_image, 2));
  surface = tomo.trws2(single(new_image), single(at_slope), single(at_weight), single(ct_slope), single(ct_weight), uint32(max_loops), uint32(bounds));

end