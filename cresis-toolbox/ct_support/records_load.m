function records = records_load(param,recs,field_names)

if ~exist('field_names','var') || isempty(field_names)
  %% Load all fields
  records_fn = ct_filename_support(param,'','records');
  records = load(records_fn);
  % Add new field "bit_mask" if missing
  if ~isfield(records,'bit_mask')
    records.bit_mask = zeros(size(records.offset));
  end
  % Remove deprecated field "surface" if present
  if isfield(records,'surface')
    records = rmfield(records,'surface');
  end
end