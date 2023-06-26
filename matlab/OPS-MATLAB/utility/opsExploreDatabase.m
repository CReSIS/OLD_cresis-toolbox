% script opsExploreDatabase
%
% Prints out the contents of the OPS database for debugging

% system: string containing accum, auth, kuband, opsuser, rds, or snow
system = 'opsuser';

query = 'SELECT * FROM pg_catalog.pg_tables';
[status,tables] = opsQuery(query);

%query = 'select schemaname, viewname from pg_catalog.pg_views where schemaname NOT IN (''pg_catalog'', ''information_schema'') order by schemaname, viewname;';
%query = 'select * from pg_catalog.pg_views';
%[status,tables] = opsQuery(query);

if 0
  %% Print all the tables
  tables_sorted = sort(tables(2,:))
  for idx = 1:size(tables,2)
    if ~strncmp(tables_sorted{idx},'pg_',3)
      fprintf('%s\n', tables_sorted{idx});
    end
  end
end

MAX_ROWS = 5;
MAX_WIDTH = 50;
for idx = 1:size(tables,2)
  if strncmp(tables{2,idx},system,length(system))
    fprintf('=============================================\n');
    fprintf('Table %s:\n', tables{2,idx})
    query = sprintf('select column_name from information_schema.columns where table_name=''%s'';', tables{2,idx});
    [status,columns] = opsQuery(query);
    
    %% Print the contents of the table out
    query = sprintf('SELECT * FROM %s LIMIT %d', tables{2,idx}, MAX_ROWS);
    [status,data] = opsQuery(query);
    
    data = data';
    
    if status == 1
      % First we convert each table entry into a string and keep track of
      % the string length
      tableContents = {};
      columnWidths = zeros(1,size(data,1));
      for row = 1:size(data,1)
        for col = 1:size(data,2)
          if isnumeric(data{row,col})
            tableContents{row,col} = sprintf('%g',data{row,col});
          elseif islogical(data{row,col})
            tableContents{row,col} = sprintf('%g',data{row,col});
          else
            tableContents{row,col} = sprintf('%s',data{row,col});
          end
          tableContents{row,col} = tableContents{row,col}(1:min(MAX_WIDTH,end));
          columnWidths(row,col) = length(tableContents{row,col});
        end
      end
      
      % We then find the longest entry in each column
      columnWidths = max(columnWidths,[],1);
      
      % We then print out the column headers with the widths being equal to the
      % longest entry
      for col = 1:length(columns)
        columnWidths(col) = max(columnWidths(col),length(columns{col}));
        fprintf(sprintf('%%-%ds',columnWidths(col)+2), columns{col});
      end
      fprintf('\n');
      
      % Next, print out the table entries
      fprintf('%s','-'*ones(1,sum(columnWidths)+2*length(columnWidths))); fprintf('\n');
      for row = 1:size(data,1)
        for col = 1:size(data,2)
          fprintf(sprintf('%%-%ds',columnWidths(col)+2), tableContents{row,col});
        end
        fprintf('\n');
      end
    end
    
    fprintf('\n');
  end
end

return
