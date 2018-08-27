function [ stats_struct ] = compute_hit_ratios( stats_struct )
%   compute the hit ratios: how much the errors (predicted value vs. ground
%   truth) are within certain numbers (e.g. 5 or 10 or 15 or ...)

%   input: 
%   stats_struct: a structure used for statistics purposes (after
%   hyperparameters tests)

    dl = length(stats_struct.total_diff);
    hit_ratios = [0 5 10 15 20 25];
    
    ratios = [];
    
    for i = hit_ratios
      if i == 0
        less_than_x = stats_struct.total_diff(abs(stats_struct.total_diff) == i);
      else
        less_than_x = stats_struct.total_diff(abs(stats_struct.total_diff) < i);
      end      
      hits = length(less_than_x);
      hit_ratio = hits/dl;
      ratios = [ratios hit_ratio];
    end    
        
    stats_struct.hit_ratios = ratios;    
end
