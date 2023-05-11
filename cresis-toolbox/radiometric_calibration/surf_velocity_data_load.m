hara;


raw_v = importdata('~/scripts/cresis-toolbox/cresis-toolbox/radiometric_calibration/GRE_G0240_2017_v.tif');

vars = fieldnames(raw_v);
for i = 1:length(vars)
    assignin('base', vars{i}, raw_v.(vars{i}));
end

