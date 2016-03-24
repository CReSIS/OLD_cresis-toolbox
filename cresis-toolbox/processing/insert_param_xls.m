% This program will not clear the entire worksheet before writing
%   If a field were deleted or elements of a nested field(wfs) were reduced
%   for a date, unwanted columns may remain on the worksheet% struct & double rows
% This program will not create a new worksheet
% If nested structure contains only 1 element, is of type
%   a_:(field)(element), and has not been marked as a special field, the
%   program will NOT find the column
% If field was 'nested structure' and is now 'double' or 'string' the
%   headers will remain, but the old data for that header and other headers
%   will be deleted.
%
%
%
% Updates param spreadsheet worksheets according to editted param structure
% in Matlab.
%   Edits rows and columns
%   Adds rows and columns
%   Retains format of unaltered param spreadsheet
%
% Inputs:
%   param_fn: string containing parameter spreadsheet filename
%   param: param structure to be imposed on parameter spreadsheet
%   worksheets: worksheet in parameter spreadsheet to be editted
%   special_fields: fields to have ar:(field)(n) format despite # of
%     elements
%
%
% Example:
%-------------------------------------------------------------------------------------------------
%   mkdir('H:\','insert_param_xls_example');
%   param_original = read_param_xls('H:\scripts\params\rds_param_2011_Greenland_P3.xls');
%   copyfile('H:\scripts\params\rds_param_2011_Greenland_P3.xls','H:\insert_param_xls_example\')
%   param_fn = 'H:\insert_param_xls_example\rds_param_2011_Greenland_P3.xls';

% 
%
%   param_edit = param_original;
%   for i = 1:numel(param_edit)
%     param_edit(i).radar.fs = i*10;
%     param_edit(i).radar.wfs(1).ref_fn(1).a = 1;
%     param_edit(i).radar.wfs(1).ref_fn(1).b = 2;
%     param_edit(i).radar.wfs(1).ref_fn(2).a = 3;
%     param_edit(i).radar.wfs(1).ref_fn(2).b = 4;
%     param_edit(i).radar.wfs(3) = param_edit(i).radar.wfs(2);
%     param_edit(i).radar.new(1).field1 = i*10;
%     param_edit(i).radar.new(1).field2 = i*20;
%   end
%   insert_param_xls(param_fn,param_edit,'radar');
%   param_new = read_param_xls(param_fn);
%-------------------------------------------------------------------------------------------------
% ispc
% Author: Jordan Sprick
%
% See also: edit_param_xls.m, insert_param_xls.m, read_param_xls.m,
% create_param_xls.m


function insert_param_xls(param_fn,param,worksheets,special_fields)
  
  [status,sheets,xlFormat] = xlsfinfo(param_fn);
  if isempty(xlFormat)
    error('Excel must be installed to use this command');
  end

  if ~ischar(param_fn)
    error('insert_param_xls:param_fn','''param_fn'' input must be a string'); 
  end
  
  if ~isstruct(param)
    error('insert_param_xls:param','''param'' input must be a struct'); 
  end

  if nargin<3 || isempty(worksheets)
    worksheets = {'cmd','vectors','records','get_heights','csarp','combine','radar'};
  elseif ischar(worksheets)
    worksheets = {worksheets};
  elseif ~iscell(worksheets)
    error('insert_param_xls:worksheet_input','''worksheet'' input must be a string or cell of strings');    
  end
  
  if nargin<4
    special_fields = [];
  elseif ischar(special_fields)
    special_fields = {special_fields};
  elseif ~(ischar(special_fields)||iscell(special_fields))
    error('insert_param_xls:special_fields_input','''special_fields'' input must be a string or cell of strings'); 
  end

  param_size = length(param);


  %% CHECKS FOR DATE_SEG REDUNDANCY

  day_segs = zeros(param_size,11);

  % loops through day_seg field of param and stores in day_segs
  for d_s = 1:param_size
    day_segs(d_s,:) = param(d_s).day_seg;
  end

  % finds unique date_segs
  [uni_rows,~] = size( unique(day_segs,'rows') );

  % checks that number of unique date_segs and total date_segs are equal
  if uni_rows<param_size
    error('insert_param_xls:redundant_day_seg','A day_seg was repeated in the param struct.');
  end



  %% READS PARAM STRUCT && WRITES TO XLS WORKSHEETS

  % loops through worksheets of spreadsheet
  for wks_idx = 1:length(worksheets)

    % This is only necessary due to the difference in 'cmd'/'command'
    worksheet = worksheets{wks_idx};    % acquires worksheet string from vector

    fprintf('Editing %s worksheet.\n',worksheet);

    % loads xls data
    try
      [ ~, ~, all_data] = xlsread(param_fn,worksheet);
    catch message
      error(message.identifier,message.message);
    end
    
    
    if any(strcmpi(worksheet,{'cmd','command'}))
      xls_version = str2double(all_data{1,2});
      if xls_version < 4
        warning('insert_param_xls:param_file_version','File version is %.1f. This function only works for 4.0 and later. Cannot write to ''command'' worksheet.');
        continue;
      elseif isnan(xls_version)
        warning('Could not find version field text in row 1, column 2 of the first worksheet.');
        continue;
      end
    end

    % finds header row (assuming 'Date' is in header column
    header_row = find(strcmp( 'Date' , all_data(1:end,1) ));

    % gets number of rows and columns for whole all_data
    col_num = size(all_data,2);

      % loops through columns of data for NaN value
      for col = 1:col_num
        data = all_data(header_row,col);
        data = data{1};
        if isnan(data)
          if col~=2
            all_data = all_data(:,1:col-1);   % deletes extra NaN columns
            break;
          end
        end
      end

    % gets number of rows and columns in trimmed all_data
    col_num = size(all_data,2);
    % creates cell to containing NEW data acquired in PARAM
    worksheet_data = [all_data(header_row:header_row+1,:) ; cell(param_size,col_num)];

    column_data = cell(param_size,1);   % cell array to contain data in column

    % DATE column
    % reads date_seg from PARAM,isolates the date,stores dates
    for r = 1:param_size
      cell_data = param(r).day_seg;
      column_data{r} = cell_data(1:8);
    end
    worksheet_data(3:end,1) = column_data;

    % SEGMENT column
    % reads date_seg from PARAM,isolates the segment,stores segments
    for r = 1:param_size
      cell_data = param(r).day_seg;
      cell_data = cell_data(10:11);
      column_data{r} = str2double(cell_data);
    end
    worksheet_data(3:end,2) = column_data;

    
    % populates data columns of worksheets
    worksheet_data = populate_worksheet_data(param,worksheet,worksheet_data,special_fields);

    % gets first cell index of header row
    xls_idx = ['A' num2str(header_row) ];

    fprintf('   Writing to %s worksheet.\n',worksheet);
    
    % writes entire worksheet array from 1st cell of header to last cell
    [~,message] = xlswrite(param_fn,worksheet_data,worksheet,xls_idx);
    % checks for writing issues
    if ~isempty(message.message)
      message.message = strrep(message.message,'\','\\');
      error(message.identifier,message.message);
    end

  end

  fprintf('Write Complete\n');
end









function worksheet_data = populate_worksheet_data(param,field,worksheet_data,special_fields,recurse)
% Extracts param structure data and stores data in existing worksheet_data
% cell array for write to xls worksheet
%
% for tx_weights in wfs of radar : does not support HANNING format
% cannot delete excess elements of structures in excel

  if nargin < 5
    recurse = 0;
  end

  headers_xls = worksheet_data(1,:);    % gets worksheet headers
  types = worksheet_data(2,:);          % gets column types
  param_size = numel(param);      % number of data columns
  special_field_flag = 0;         % flag for special_fields
  
  % Finds # of elem of struct to loop through and sizes of elements of
  %   struct
  elements = 0;   %initializes elements
  size_array = zeros(1,param_size);   %stores sizes of each element
  for p = 1:param_size
    elements_temp = numel(param(p).(field));  % number of elements of data
    size_array(p) = elements_temp;
    if elements_temp>elements;
      elements = elements_temp;     % number of elements of struct to loop through
    end
  end
  
  % loops through each element of struct str(1),str(2)...
  for e = 1:elements
    
    clear n_param
    %%  Brings out nested elements to form new structure for easier access
    % creates empty structure of same format
    fns = cell(0);
    for p = 1:param_size
      if ~isempty(param(p).(field))
        try
          fn = fieldnames(param(p).(field)(e));
          fn_mem_data = ismember(fn,fns);
          fns_new = fn(~fn_mem_data);
          if ~isempty(fns_new)
            for i = 1:length(fns_new)
              % fills preceding elements with fields of empty data
              n_param(param_size+1).(fns_new{i}) = 1;
            end
          end
          fns = {fns{:} fns_new{:}};
        catch
        end
      end
    end
    n_param(param_size+1) = [];   % erases dummy element
    
    % fills elements of structure with data
    for p = 1:param_size
       if ~isempty(param(p).(field))
         try
          fn = fieldnames(param(p).(field)(e));
          for i = 1:length(fn)
            n_param(p).(fn{i}) = param(p).(field)(e).(fn{i});
          end
         catch
         end
       end
    end
    %%
    
    fields_nested = fieldnames(n_param(1));
    column_data = cell(param_size,1);    % creates empty cell to store row data
    
    % checks if field should be in ar:_ format
    if ~isempty(find(strncmp(field,special_fields,length(field)))) && recurse
      special_field_flag = 1;
    end
    
    % JORDAN: eventually change to find actual type for NEW columns
    % Should be fine for adding additional elements of structure
    type = 'r';


    %% Handles ar:_ formatted columns
    if elements > 1 || special_field_flag
      % ar:struct(e)
      type = sprintf('ar:%s(%d)',field,e);
      % finds columns with matching ar:_(e) type
      type_matches = find(strncmpi(type,types,length(type)));
      % gets headers of columns
      headers = headers_xls(type_matches);
      
      % deal with SIZE column first
      % SIZE column must be first so read_param_xls works correctly
      if e == 1
        % ar:struct
        size_type = sprintf('ar:%s',field);
        % finds column with matching size_type
        size_type_matches = find(strncmpi(size_type,types,length(size_type)+1));
        % if SIZE column exists in xls worksheet
        if ~isempty(size_type_matches)
          % add to type_matches and headers
          type_matches(end+1) = size_type_matches;
          headers{end+1} = 'size';
        end
        % adds fields so that SIZE column will be stored in next loop
        fields_nested = {'size',fields_nested{:}};
        % stores sizes of 
        for p = 1:param_size
          n_param(p).size = size_array(p);
        end
      end
    end
    %% 
    % Loops through fields of structure
    for f = 1:length(fields_nested)
      
      field_nested = fields_nested{f};
      
      % Gets class of all field data
      type_samples = cell(0);
      for r = 1:param_size
        if ~isempty(n_param(r).(field_nested))
          type_samples{end+1} = class(n_param(r).(field_nested));
        end
      end
      % Gets unique classes of field data
      if isempty(type_samples)
        type_samples_unique = 'double';
      else
        type_samples_unique = unique(type_samples);
      end
      
      %% Checks for existing column headers
      if ~recurse
        % checks for xls columns with matching header
        column = find(strcmp(headers_xls,field_nested));
        % checks for columns of 'ar:_' type with same header
        type_matches = find(strncmpi('a',types,1));
        % finds columns NOT of 'ar:_' type
        column = column(~ismember(column,type_matches));

      elseif elements > 1 || special_field_flag
        % checks for xls columns with 'ar:_' headers
        header_idx = find(strcmp(headers,field_nested));
        % types_matches has already been acquired while looking for
        %   existing 'size' column
        column = type_matches(header_idx);
        if isempty(column)
          header = strcat(field,'.',field_nested);
          column = find(strcmp(headers_xls,header));
            if ~isempty(column)
              error('insert_param_xls:header_format',sprintf('''%s'' column in worksheet has wrong format for field of multiple elements.',field));
            end
        end

      elseif elements==1
        % checks for xls columns with 'struct.field' headers
        header = strcat(field,'.',field_nested);
        column = find(strcmp(headers_xls,header));
        if isempty(column)
          type = sprintf('ar:%s',field);
          type_matches = find(strncmpi(type,types,length(type)));
          if ~isempty(type_matches)
            error('insert_param_xls:header_format','Field column in worksheet has wrong format for field of 1 elements that isn''t specified as a special field.');
          end
        end
      end
      
      if length(column)>1
        error('insert_param_xls:header_redundancy','Header was found multiple times in worksheet');
      end
      %%
      
      %% Recurses or creates a new column
      % Found no existing column header in xls worksheet
      if isempty(column)
        
        % Makes sure all field data is of struct class
        % Makes sure program hasn't already recursed
        if length(type_samples_unique) == 1 && strcmp(type_samples_unique,'struct') && ~recurse
          % RECURSES to inside of data inside of doubly nested structure
          worksheet_data = populate_worksheet_data(n_param,field_nested,worksheet_data,special_fields,1);
          continue;
        end
        
        % Create new column to contain field data
        % Finds next empty column
        column = size(worksheet_data,2) + 1;
        % Stores field type for new column
        if strcmp(field_nested,'size')
          worksheet_data(1:2,column) = {field_nested,size_type};
        elseif ~recurse || elements > 1 || special_field_flag
          worksheet_data(1:2,column) = {field_nested,type};
        elseif recurse
          worksheet_data(1:2,column) = {strcat(field,'.',field_nested),type};
        end
        types{end+1} = type;
      end
      %%

      % loops though elements(rows) of field(column)
      for r = 1:param_size;

        if ~strncmp('struct',class(n_param(r).(field_nested)),6)
          % Stores non-struct data raw
          cell_data = n_param(r).(field_nested);
        else
          % Re-formats struct data to string
          cell_data = struct_to_matlab_cmds_edit(n_param(r).(field_nested),1);
        end

              cell_num = 1;
              % Checks if data is of cell class (requires special
              %   attention)
              if iscell(cell_data)
                % # elements in cell array
                cell_num = numel(cell_data);
              end

              % - runs once for non-cell data
              % - iterates for every element for cell array data
              for cell_idx = 1:cell_num

                if iscell(cell_data)
                  % Gets indexed element of cell array
                  original_data = cell_data{cell_idx};
                  % Begins string containing cell array data
                  if cell_idx == 1;
                    final_data = '{';
                  end
                else
                  original_data = cell_data;
                end

                % Reformats cell data need special formatting
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Empty COMMAND worksheet cell (not necessary but looks better)
                if strcmp(field,'cmd') && ~isempty(original_data) && islogical(original_data)
                  if original_data == 0
                    formatted_data = [];
                  % fixes xls 1->'true'
                  elseif original_data == 1;
                    formatted_data = 1;
                  end

                % binary
                elseif strcmpi(types{column},'b')
                  if original_data == 1;
                    formatted_data = 1;
                  else
                    formatted_data = 0;
                  end

                %adc_gains
                elseif strcmpi(field_nested,'adc_gains') && ~isempty(original_data)
                  x = -(20*log10(original_data)-72);
                  x_unique = unique(x);
                  x_num = length(x_unique);
                  if x_num>1
                    mat_str = mat2str(x);
                    formatted_data = sprintf('10.^((72-%s)/20);', mat_str);
                  else
                    formatted_data = sprintf('10.^((72-%d*ones(1,%d))/20);', x_unique, length(original_data));
                  end

                % Matlab matrix format
                elseif isnumeric(original_data) && numel(original_data)>1
                  formatted_data = mat2str(original_data);

                % Matlab function_handle format
                elseif isa(original_data,'function_handle')
                  formatted_data = func2str(original_data);
                  formatted_data = ['@' formatted_data];

                % Matlab inline format
                elseif isa(original_data,'inline')
                  formatted_data = char(original_data);
                  formatted_data = ['inline(''' formatted_data ''')'];

                % Matlab String format (not 't' type)
                %   makes sure string is not re-formatted structure
                %   string
                elseif ischar(original_data) && ~isempty(original_data) && ~strncmp('struct',original_data,6)...
                       && ~strcmpi(field,'cmd') && (strcmp(types{column},'r') || strncmp('ar',types{column},2))
                  try
                    eval([original_data ';']);
                  catch
                    % adds quotes if string cannot be eval()
                    formatted_data = sprintf('''''%s''',original_data);
                  end

                % XLS cell is INF
                elseif isinf(original_data)
                  formatted_data = 'Inf';

                % logical of 'r' type
                elseif ~isempty(original_data) && islogical(original_data) && strcmp(types{column},'r')
                  if original_data == 1
                    formatted_data = '''true';
                  else
                    formatted_data = '''false';
                  end

                % no re-formatting
                else
                  formatted_data = original_data;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % forms string for cell class data
                if iscell(cell_data)
                  if isnumeric(formatted_data)
                    formatted_data = num2str(formatted_data);
                  end
                  final_data = sprintf('%s,%s',final_data,formatted_data);
                  if cell_idx == cell_num
                    final_data = sprintf('%s}',final_data([1,3:end]));
                  end

                else
                  final_data = formatted_data;
                end

              end

        % writes data cell to column vector
        column_data{r} = final_data;
      end

      % Removes the created 'size' field to correctly read other
      % following elements
      if strcmp(field_nested,'size')
        n_param = rmfield(n_param,'size');
      end

      % Writes column vector to worksheet data cell array
      worksheet_data(3:end,column) = column_data;
    end
  end
end