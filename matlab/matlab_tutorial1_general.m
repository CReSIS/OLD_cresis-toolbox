% =========================================================================
%                      CReSIS MATLAB TUTORIAL
% =========================================================================
%
% 1. Follow steps here:
%    https://wiki.cresis.ku.edu/cresis/MATLAB_Tutorial
% 2. Run this tutorial by pressing F5 (choose "Change Folder" if asked)
%    Follow along in the terminal (command window) and editor
%
% Authors: Audrey Evans, Megan Metz, John Paden, Kyle Purdon, Levi Sedlock,
%  Jordan Sprick

clc; clear; close all;
format compact;
fprintf('Available tutorials:\n');
fprintf('1: Variables\n');
fprintf('2: Structures and Cell Arrays\n');
fprintf('3: Conditionals and loops\n');
fprintf('4: file i/o\n');
fprintf('5: visualization\n');
fprintf('6: mapping toolbox\n');
fprintf('7: vectorization\n');
fprintf('8: debugging\n');
fprintf('9: signal processing\n');
fprintf('10: signal processing 2\n');
fprintf('11: signal processing 3: s-parameters\n');
fprintf('12: Layer Plots\n');

done = false;
while ~done
  try
    section_number = input('Enter tutorial to run [1-12]: ');
    section_number = round(section_number(1));
    if section_number >= 1 && section_number <= 12
      done = true;
    end
  end
end

if any(section_number == [4 5 6 11 12])
  global gdata_folder_path;
  if isempty(gdata_folder_path)
    gdata_folder_path = 'C:\tmp\MATLAB_tutorial_files\';
  end
  fprintf('\nDefault Data Path: %s\n', gdata_folder_path);
  fprintf('Where did you place the data? Press enter to accept the default.\n');
  done = false;
  while ~done
    try
      data_folder_path = input(': ','s');
      if isempty(data_folder_path)
        data_folder_path = gdata_folder_path;
      end
      if ~exist(data_folder_path,'dir')
        fprintf('Does not exist: %s\n', data_folder_path);
      else
        gdata_folder_path = data_folder_path;
        done = true;
      end
    end
  end
end

clc; clear done;
edit MATLAB_Tutorial.m
fprintf('Entering debug mode. Use the Matlab editor to view the code as it runs.\n');
fprintf('The green arrow in the editor indicates the line that will be run next.\n')
fprintf('In debug mode, you can step through the code line by line.\n');
fprintf('You can also enter any ordinary Matlab commands.\n');
fprintf('Debug Commands:\n');
fprintf('  STEP FORWARD ONE LINE:\n');
fprintf('    Press F10, run "dbstep", or menu item "Editor-->Step"\n');
fprintf('  RUN THE REST OF THE FUNCTION WITHOUT STOPPING:\n');
fprintf('    Press F5, run "dbcont", or menu item "Editor-->continue"\n');
fprintf('  QUIT:\n');
fprintf('    Press Shift + F5, run "dbquit", or menu item "Editor-->quit debuging"\n');

%% Section 1: variables
if section_number == 1
  fprintf('\n=========================================================\n');
  fprintf(' Section 1: Variables\n');
  fprintf('=========================================================\n');
  fprintf('Types of variables tutorial:\n');
  fprintf('  http://www.mathworks.com/help/toolbox/eml/ug/bq2l_yi.html\n');
  fprintf('Overview:\n');
  fprintf('  http://www.mathworks.com/videos/introducing-matlab-fundamental-classes-data-types-68991.html\n');
  keyboard
  
  % TYPES OF VARIABLES
  % =======================================================================
  % http://www.mathworks.com/help/matlab/data-types_data-types.html
  
  % Variable assignments
  a = 3
  b = 2

  
  % OPERATIONS ON VARIABLES
  % =======================================================================
  % http://www.mathworks.com/help/matlab/operators-and-elementary-operations.html
  % Using variables (all the regular mathematical operations are available)
  c = a+b
  c = a-b
  c = a*b
  c = a/b
  c = a^b
  

  % MATRICES
  % =======================================================================
  % All variables are treated as matrices in Matlab.
  % - A scalar (single number) is treated as a 1 by 1 matrix.
  % - Matrices can have any number of dimensions
  % - The first dimension is the row, the second dimension is the column
  % - New rows are started with ";"
  
  a = [1 2 3 4]
  b = [1; 2; 3; 4]
  
  
  % MATRIX OPERATIONS
  % =======================================================================
  % Matrix operations
  % - Multiply each element of a by 2:
  c = a*2
  % - Divide each element of b by 2:
  c = b/2
  % - Element-wise addition and subtraction
  c = a+[2 2 2 2] % Add 2 to each element
  c = a+[0 0 0 1] % Add 1 to the last element
  % - Element-wise multiply
  c = a.*a
  % - Element-wise divide
  c = a./a
  % - Matrix multiply (linear algebra)
  c = a*b
  % - Matrix transpose (rows are turned into columns and vice versa)
  c = a.'
  % - Logical operations
  c = a < 2
  c = a <= 2
  c = a == 2
  c = a > 2
  c = a >= 2
  c = ~(a >= 2) % ~ is "Not"
  
  % Some common matrices can be created easily
  a = 1:4
  a = 0:0.5:4
  a = linspace(0,10,5)
  
  
  % INDEXING MATRICES 1
  % =======================================================================
  % Indexing: How to access individual elements of the array
  e = [1 2 3; 4 5 6]
  
  e(1,1)
  e(2,1)
  e(1,3)
  e(2,2)
  % Shortcut to the last element in each dimension
  e(1,end) % Same as e(1,3)
  e(end,1) % Same as e(2,1)
  e(end,end) % Same as e(2,3)
  % How many rows
  size(e,1)
  % How many columns
  size(e,2)
  % How many elements
  numel(e)
  
  
  % INDEXING MATRICES 2
  % =======================================================================
  % Indexing: How to access many elements of the array
  % First column:
  e(:,1)
  % First row:
  e(1,:)
  % First two elements of the first row (use another matrix to indicate
  % which elements)
  e(1,[1 2])
  % Using another matrix variable to index the matrix
  idxs = [1 2];
  e(:,idxs)
  % Last two elements of the first row (use another matrix to indicate
  % which elements)
  e(1,[2 3])
  e(1,[end-1 end])
  % Combine two matrices: concatenate in column dimension
  [e e]
  % Combine two matrices: concatenate in row dimension
  [e; e]
  % Concatenate a subset of the matrix
  [e(:,1) e(:,1)]
  
  
  % INDEXING ELEMENTS 3
  % =======================================================================
  % Adding elements to a matrix by assigning to indices beyond the existing
  % elements of the matrix
  e(3,:) = [7 8 9]
  e(end+1,:) = [10 11 12]
  
  
  % SUPPRESS OUTPUT
  % =======================================================================
  % Suppress the printout of the variable contents with a semi-colon at the
  % end of the line
  e
  e; % Notice that nothing gets printed in the terminal
  2+2
  2+2; % Notice that nothing gets printed in the terminal
  100-24; % Notice that nothing gets printed in the terminal

  
  % EMPTY VARIABLES
  % =======================================================================
  % Empty variable
  e = []
  size(e,1)
  size(e,2)
  numel(e)

  
  % STRING or TEXT
  % =======================================================================
  % Strings: variables which contain text (these are matrices of ascii
  % charcters in Matlab)
  
  % Use single quotes to start and stop a string
  s1 = 'Hello'
  s2 = 'there'
  s3 = 'to put a quote in a string, use ''two quotes'' in a row'

  % Combine strings
  s = [s1 s2]
  s = [s1 ' ' s2]

  % Indexing works the same
  s1(1,[1 2 3 4])
  s1(1,1:4)
  s1(1,[end-3 end-2 end-1 end])
  s1(1,end-3:end)

  
  % CLASS OF A VARIABLE
  % =======================================================================
  % class: The type of the variable
  a = 3.4235 % double: which means a regular number
  class(a)
  a = true % logical: true or false
  class(a)
  a = 'text' % string (known as "char" in Matlab which is short for "ASCII character")
  class(a)

  
  % LOGICAL INDEXING
  % =======================================================================
  % "Logical indexing" is a way to get elements of a matrix which satisfy
  % certain properties.
  numbers = [3 23 4 3 55 2 1 4]
  
  % Two ways to grab the first 3 numbers:
  numbers([1 2 3]) % Regular way
  numbers(([true true true false false false false false])) % Logical indexing
  
  % Grab all numbers less than 10
  numbers(numbers < 10)
  % Set all numbers greater than 10 to 10
  numbers(numbers > 10) = 10

  % Find the index for all elements that are less than 10
  find(numbers < 10) % "find" returns the indices of all "true" elements
  
  fprintf('Section 1 Complete.\n');
  return
end

%% Section 2: Structures and Cells
if section_number == 2
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 2: Structures\n');
  fprintf('=========================================================\n');
  fprintf('\n')
  fprintf(' Please watch this video to introduce structures and cells:\n')
  fprintf(' http://www.mathworks.com/videos/introducing-structures-and-cell-arrays-68992.html\n')
  fprintf(' follows are two example of how to create/access structures and cells\n')
  fprintf('\n')
  keyboard
  

  % STRUCTURE
  % =======================================================================
  % Structures are another class that let you organize variables together.
  % Each "sub-variable" is called a "field" in the structure. The
  % sub-variables can be any class (including another structure)
  
  % Create a "student" structure with four fields
  clear student % Erase the variable first
  student.first_name = 'Kyle'
  student.last_name = 'Purdon'
  student.email = 'kylepurdon@mail.com'
  student.age = 20
  % Add a second "student" to the structure
  student(2).first_name = 'Trey';
  student(2).last_name = 'Stafford';
  student(2).email = 'treystafford@mail.com';
  student(2).age = 19;
  % Add a third "student" to the structure (but leave out age)
  student(3).first_name = 'First';
  student(3).last_name = 'Last';
  student(3).email = 'pirate@mail.com';
  % Change the name of the last structure
  student(end).first_name = 'Davy';
  student(end).last_name = 'Jones';
  
  student
  student(2).first_name
  student(3).first_name
  student(1).age
  student(3).age % Undefined fields will be empty until defined
  
  
  % CELLS
  % =======================================================================
  % Cells are another class that let you organize variables, but each
  % element is indexed by a number rather than a name
  mycell = {student [1 2; 3 4] 'Cell Arrays are Awesome!'}

  % Accessing a single element or elements of a cell array will return
  % another cell array containing the indexed elements.
  mycell(1)
  mycell(2)
  mycell(3)
  mycell(1:2)
  
  % To access the CONTENTS of the cell array use {} instead of ()
  mycell{1}
  mycell{2}
  mycell{3}
  mycell{1}(2).last_name % Last name of the second student
  
  
  % Cell arrays are most common when making lists of strings
  filenames = {'data_file1.txt','data_file2.txt','data_file23.txt'}
  filenames{1}
  filenames{3}
  % Add a new filename to the list
  filenames{end+1} = 'new_datafile.txt'
  
  fprintf('Section 2 Complete.\n');
  return
end

%% Section 3: conditionals and loops
if section_number == 3
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 3: Conditionals and Loops\n');
  fprintf('=========================================================\n');
  
  keyboard

  % IF
  % =======================================================================
  % If statements test a condition and only run the commands inside if it
  % is true.
  number = 12;
  if number > 10
    fprintf('The number is greater than 10.\n');
  end

  % Else statements can be used to run commands with the condition fails
  number = 5;
  if number > 10
    fprintf('The number is greater than 10.\n');
  else
    fprintf('The number is less than 10.\n');
  end
  
  % Elseif statements allow multiple conditions to be tested
  if number < 5
    fprintf('The number is less than 5.\n');
  elseif number == 5
    fprintf('The number is 5.\n');
  else
    fprintf('The number is greater than 5.\n');
  end
  
  % "any" and "all" are commonly used functions with IF and allow you to
  % test many elements in a matrix at once to see if any or all are true.
  if any([false false true])
    fprintf('At least one was true.\n');
  end
  
  if all([false false true])
    fprintf('All were true.\n');
  else
    fprintf('At least one was not true.\n');
  end
  
  if all([true true true])
    fprintf('All were true.\n');
  end
  
  % If you need to combine conditions use "&&" for and "||" for or
  condition1 = true;
  condition2 = false;
  if condition1 && condition2
    fprintf('Both conditions are true.\n');
  elseif condition1 || condition2
    fprintf('At least one condition is true.\n');
  end

  
  % SWITCH
  % =======================================================================
  % Switch statements are just like if...elseif...else statements. They are
  % more limited in how they can be used, but are more organized.
  
  stop_light = 'yellow'
  switch stop_light
    
    case 'red'
      fprintf('The light is red');
      
    case 'yellow'
      fprintf('The light is yellow');
      
    case 'green'
      fprintf('The light is green');
      
    otherwise
      fprintf('The light is not working...');
  end
  
  
  % FOR LOOPS
  % =======================================================================
  % For loops are used to execute the same piece of code many times for
  % different variable values.  The statements of code in a for loop
  % execute for each column of the for loop variable. Each time it executes
  % the variable will be set to the column corresponding to the number of
  % the times the for loop has executed
  
  % Print out the first four squares
  for x = 1:4
    x*x
  end

  % Print out columns of matrix
  e = [1 2 3; 4 5 6]
  for column = e
    column
  end

  % Print out a list of strings in a cell array
  colors = {'red','blue','green'}
  for color = colors
    color
  end

  % What if the cell array is a single column?
  colors = {'red'; 'blue'; 'green'}
  for color = colors
    color
  end
  % Using indices to access the individual elements
  for index = 1:numel(colors)
    index
    color = colors{index}
  end

  % Compute 4-factorial (i.e. 1*2*3*4)
  number = 4;
  result = 1;
  for n = 1:number
    result = result*n;
  end
  fprintf('The factorial of %d is %d.\n',number,result)
  
  % What is the above doing? The loop will run from 1 to "number". If the
  % number is 4 the values each run of the loop are shown below.
  %
  % n=1  fact = 1*1 = 1
  % n=2  fact = 1*2 = 2
  % n=3  fact = 2*3 = 6
  % n=4  fact = 6*4 = 24

  
  % WHILE LOOPS
  % =======================================================================
  % The statements in a while loop run until the while-loop condition
  % fails.
  
  % Another way to computer factorial
  number = 4;
  result = 1;
  while number > 1
    result = result * number
    number = number - 1
  end

  
  % BREAK AND CONTINUE
  % =======================================================================
  % break: Stop executing in the middle of a loop
  % continue: Skip to the next iteration
  for idx = 1:10
    if idx < 3
      continue;
    end
    idx
    if idx > 5
      break;
    end
  end
  
  
  fprintf('Section 3 Complete.\n');
  return
end

%% Section 4: File Input and Output
if section_number == 4
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 4: File IO\n');
  fprintf('=========================================================\n');
  
  fprintf(' Watch this video for an introduction to file I/O\n')
  fprintf(' http://www.mathworks.com/videos/importing-data-from-files-68988.html\n')
  fprintf('\n');
  fprintf(' Below we will give 2 file reading examples.\n')
  fprintf(' 1) Read a standard CReSIS CSV\n')
  fprintf(' 2) Read a CReSIS MAT file\n')
  fprintf('\n')
  
  keyboard
  
  % READ CRESIS CSV FILE (CONTAINS ICE SURFACE AND ICE BOTTOM)
  % =======================================================================
  
  % file_path: string containing the path to the file
  % fullfile: Matlab's operating system independent file path creation.
  %   This function is the best way to build paths in Matlab.
  file_path = fullfile(data_folder_path,'CRESIS_CSV_DATA.csv')
  
  % Open the file in the editor to view the contents
  %  - This file is a "comma separated variable" or csv file
  %  - The field are: LAT,LON,UTCTIMESOD,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY
  edit(file_path)
  
  % Open a file identifier to the file
  fid = fopen(file_path);
  
  % Scan and load the data from file using a specified format
  % - The format string tells how the file is formatted. For this
  %   case it is '%f%f%f%f%f%f%f%f%s'.
  %   o %f is for a number
  %   o %s is for a string
  %   o this format is 8 numbers followed by a string on each line of the
  %     file
  % - The other arguments tell textscan other file format information
  %   that it needs to read in the file.
  %   o delimiter: there are commas between each field
  %   o headerlines: there is 1 line that it is header and should be
  %     ignored
  csv_data = textscan(fid,'%f%f%f%f%f%f%f%f%s','delimiter',',','headerlines',1);
  
  % Close the data pointer (DONT FORGET THIS)
  fclose(fid);
  
  % Note that csv_data is a cell array
  % csv_data{1} is the first column of data from the file, csv_data{2} is
  % the second column of data, and so on.
  csv_data
  
  % Using multiple output arguments with a cell array:
  % When multiple cell contents are accessed at once, each represents a
  % different output. These outputs are assigned by placing a list of
  % output variable names in square brackets separated by commas like this:
  [LAT,LON,UTCTIMESOD,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY] = csv_data{:};
  
  % This file contains ice thickness data from a flight track
  h_fig = figure; % Create a new figure for plotting
  plot(LON, LAT); % plot the longitude and latitude
  xlabel('Longitude (deg)'); % add an xlabel to the plot
  ylabel('Latitude (deg)'); % add a ylabel to the plot
  
  clf(h_fig);
  scatter(LON, LAT, [], THICK); % same plot as before, but color is ice thickness
  xlabel('Longitude (deg)'); % add an xlabel to the plot
  ylabel('Latitude (deg)'); % add a ylabel to the plot
  h_colorbar = colorbar;
  set(get(h_colorbar,'YLabel'),'String','Ice thickness (m)');
  
  clear fid csv_data h_colorbar h_fig LAT LON UTCTIMESOD THICK ELEVATION FRAME SURFACE BOTTOM QUALITY;
  
  
  % READ CRESIS LAYER DATA FILE (CONTAINS ICE SURFACE AND ICE BOTTOM)
  % =======================================================================
  % The format is explained here:
  % https://wiki.cresis.ku.edu/cresis/Layer_File_Guide
  
  % Set the file path
  file_path = fullfile(data_folder_path,'CRESIS_MAT_DATA.mat')
  
  % To demonstrate the file loading, we will print out our current variables
  whos
  
  % Load the Matlab .mat file (Matlab's proprietary format) directly into
  % your current workspace. This method of loading will load all the
  % variables stored in the file and create variables with those same names
  % (potentially overwriting variables you already have defined).
  load(file_path)
  whos
  
  clear Elevation GPS_time Latitude Longitude layerData;
  
  % The previous method tends to clutter the workspace, so we generally use
  % the following to load the variables into a structure
  mat_data = load(file_path)

  % To load a few specific variables
  mat_data = load(file_path,'Latitude','Longitude')

  % "whos" can be use in function form:
  mat_data_info = whos('mat_data')
  
  % This file contains ice thickness data from a flight track
  h_fig = figure; % Create a new figure for plotting
  plot(mat_data.Longitude, mat_data.Latitude); % plot the longitude and latitude
  xlabel('Longitude (deg)'); % add an xlabel to the plot
  ylabel('Latitude (deg)'); % add a ylabel to the plot  
  
  
  % SAVE
  % =======================================================================
  % Save is a very simple command that saves variables and
  % workspaces in mat or ascii format. See "help save" for more information.
  file_path = fullfile(data_folder_path,'CRESIS_MAT_DATA_OUT.mat')
  
  % Saves the struct variable into a file
  save(file_path,'mat_data')
  
  % Saves the fields of the struct variable into a file
  save(file_path,'-struct','mat_data')
  
  % Saves on the specified fields of the struct variable into a file
  save(file_path,'-struct','mat_data','Latitude','Longitude')
  
  % Matlab's new HDF5 format is version 7.3. This has the best format for
  % sharing binary data since HDF5 is a standard format read by many
  % libraries.
  save(file_path,'-v7.3','-struct','mat_data')

  % Saves the entire workspace with all of the variables... sometimes
  % useful for debugging.
  save(file_path)

  
  % ANOTHER CSV EXAMPLE
  % =======================================================================
  % The file looks like this:
  %
  % NAME,USERID,GRADE
  % Joe,9551,F
  % Bob,9552,A
  % Jane,9553,A
  % Mike,9554,B
  % Sally,9555,C
  % Craig,9556,B
  
  % Set the file path
  file_path = fullfile(data_folder_path,'GRADES_TXT_DATA.txt');
  
  fid = fopen(file_path);
  grade_data = textscan(fid,'%s%d%s','Headerlines',1,'Delimiter',',');
  fclose(fid);
  
  [NAME,USERID,GRADE] = grade_data{:};
  
  % Get a list of just the unique grades:
  unique_grades = unique(GRADE)

  % Remove identical names in a second list
  names_only_in_first_list = setdiff(NAME(1:4),NAME(2:5))

  % Find names in both lists
  names_in_both_lists = intersect(NAME(1:4),NAME(2:5))
  
  % Next lets get the names of the students that recieved "A" using logical
  % indexing. There are several steps here:
  % - Inline function: @(grade) strcmpi(grade,'A')
  %   An "inline" function is one that is defined in the middle of the
  %   source code. This function takes one argument and compares that to
  %   the string 'A'.
  % - strcmpi: compares two strings and ignores case. True is returned
  %   when they are equal and false when not equal.
  % - cellfun: Applies our inline function to every cell entry in GRADE and
  %   returns a matrix containing the outputs of each of these comparisons.
  %   grade_mask will be the same size as GRADE.
  % - NOTE: "grade" and "GRADE" are different variables
  grade_mask = cellfun(@(grade) strcmpi(grade,'A'),GRADE)
  
  % Some functions will handle cell array inputs automatically. strcmpi
  % allows a cell array to be passed into the first argument and it will
  % then do the same operation as the above cellfun. Note: not all
  % functions support this interface and cellfun is required.
  grade_mask = strcmpi(GRADE,'A')
  
  % A "logical" array is 1's and 0's. A "1" represents a true condition a "0"
  % a false condition. The above will give us an array of 1's and 0's the
  % length of the grades array with 1's at the indexes of the A's. We can use
  % this index to get the name/s of the students with "A" for a grade.
  students = NAME(grade_mask)      % Keep only names that have "1" in the "grade_mask"
  
  % Create a new CSV file with just these students
  file_path = fullfile(data_folder_path,'GRADES_TXT_DATA_A.txt');
  
  fid = fopen(file_path,'w+');   % The 'w+' means fopen will create a new file or discard existing contents of a file.
  
  % Print header
  fprintf(fid,'%s,%s,%s\n','NAME','USERID','GRADE');
  
  % Ensuring a vector is a row vector and not a column vector:
  grade_mask
  grade_mask(:) % Always a column vector
  grade_mask(:).' % Transpose of a column vector is always a row vector
  
  % Now we can use a loop to print the data one line at a time.
  for student_idx = find(grade_mask(:).')
    fprintf(fid,'%s,%d,%s\n',NAME{student_idx},USERID(student_idx),GRADE{student_idx});
  end
  
  % Close the file
  fclose(fid);
  
  fprintf('Section 4 Complete.\n');
  return
end

%% Section 5: Visualization
if section_number == 5
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 5: Visualization\n');
  fprintf('========================================================\n')
  
  fprintf(' Here are a few videos to watch before this section:\n')
  fprintf(' http://www.mathworks.com/videos/using-basic-plotting-functions-69018.?s_tihtmld=srchtitle\n')
  fprintf(' http://www.mathworks.com/videos/visualizing-data-with-matlab-68917.html\n')
  fprintf('\n');
  fprintf(' Figures and Handles\n')
  fprintf(' ---------------------------------------------------------\n')
  fprintf('\n');
  fprintf(' Note the figure and handle are the basic manipulators for plots and\n')
  fprintf(' figures. Read through the documentation below for detailed\n')
  fprintf(' information...\n')
  fprintf('\n');
  fprintf(' http://www.mathworks.com/help/techdoc/learn_matlab/f3-15974.html\n')
  fprintf('\n');
  fprintf(' Let''s create a few figures...\n')
  keyboard

  % CREATING A FIGURE
  % =======================================================================
  figure(1) % Creates a new figure with the figure number 1
  close(1) % Close the figure
  
  % Creates a figure with a custom title and stores the figure in handle h.
  h = figure('Name','MyPlot','NumberTitle','off');
  delete(h); % Closes the figure without calling the close function handle
  
  
  % PLOT AND STEM
  % =======================================================================
  time = 0:.01:2;
  f = cos(2*pi*2*time);
  
  figure(1);      % creates a new figure with figure number 1
  plot(time,f);   % plot a cosine curve on figure 1
  
  figure(2);      % make another figure
  stem(time,f);   % stem a cosine curve on figure 1

  figure(1);      % since figure number 1 exists already, this just selects the figure
  % Add another plot to the figure. Make this plot red.
  % - See help plot for all the variations of color and symbol.
  hold on;
  plot(time,-f,'r');   % plot a cosine curve on figure 1
  hold off;
  
  close([1 2]);
  
  % ACCESSING HANDLES
  % =======================================================================
  time = 0:.01:2;
  f = cos(2*pi*2*time);
  
  h_fig = figure; % Create a new figure
  h_axes = axes; % Create a new axes
  h_plot = plot(time,f,'Parent',h_axes); % plot a cosine curve in the new axes
  get(h_plot) % List all the plot handle fields
  set(h_plot) % List all the plot handle fields and the options allowed
  set(h_plot,'YData',-f) % Overwrite the previous plot data
  set(h_plot,'LineWidth',2,'Marker','o','MarkerSize',12) % Overwrite the previous plot data
  hold(h_axes,'on');
  h_plot = plot(time,f,'Parent',h_axes,'Color','red');
  h_plot = plot(time,1.3*f,'Parent',h_axes,'Color',[0 1 0]);
  delete(h_fig);
  

  % PLOT CRESIS LAYER DATA ONTO AN ECHOGRAM
  % =======================================================================
  % Plot a CReSIS elevation on a time axis
  file_path = fullfile(data_folder_path,'CRESIS_MAT_DATA.mat');
  layer_data = load(file_path);
  file_path = fullfile(data_folder_path,'CRESIS_MAT_ECHODATA.mat');
  echo_data = load(file_path);
  
  % Plot the image using a base 10 log plot db()
  figure;
  imagesc(db(echo_data.Data,'power'));
  xlabel('Range lines');
  ylabel('Range bin');
  
  % Remove bad data at the bottom of the image
  echo_data.Data = echo_data.Data(1:1200,:);
  echo_data.Time = echo_data.Time(1:1200);
  
  % Plot using default x-axis (column number) and two way travel time for
  % the y-axis
  clf;
  imagesc([],echo_data.Time*1e6,db(echo_data.Data,'power'))
  hcolor = colorbar;
  set(get(hcolor,'YLabel'),'String','Relative power (dB)');
  xlabel('Range lines');
  ylabel('Two way travel time (us)');
  
  % Change the colormap
  colormap(gray); % 64 color gray scale (64 is the default)
  colormap(gray(256)); % 256 color gray scale
  colormap(1-gray(256)); % invert gray scale
  colormap(hsv(256)); % hue only colormap (used to plot phase/angle)
  colormap(jet(256)); % "jet" colormap (red is big, blue is small)
  
  % Add the layer data
  hold on;
  plot(layer_data.layerData{1}.value{1}.data*1e6,'k')
  plot(layer_data.layerData{1}.value{2}.data*1e6,'k')
  plot(layer_data.layerData{2}.value{1}.data*1e6,'k')
  plot(layer_data.layerData{2}.value{2}.data*1e6,'k')
  
  % Get the current axes
  h_axes = gca;
  
  % Create an identical figure without layers
  figure;
  imagesc([],echo_data.Time*1e6,db(echo_data.Data,'power'))
  hcolor = colorbar;
  set(get(hcolor,'YLabel'),'String','Relative power (dB)');
  xlabel('Range lines');
  ylabel('Two way travel time (us)');
  colormap(1-gray(256)); % invert gray scale
  
  % Get the new axes
  h_axes(end+1) = gca;
  
  % Link the zoom for the two axes
  linkaxes(h_axes,'xy');
  
  % Before closing: try zooming in one or the other figure and see how they
  % are linked together.
  close all;
  
  
  % SUBPLOTS
  % =======================================================================
  
  % Sub plots allow us to use a single figure window to plot multiple
  % objects on different axes.
  
  file_path = fullfile(data_folder_path,'CRESIS_MAT_ECHODATA.mat');
  echo_data = load(file_path);
  
  % What is we want to plot the echogram, surface, bottom, and time domain
  % on the same plot, but seperated?
  
  %                       |----------|----------|
  %                       |          |          |
  %                       |   p=1    |   p=2    |
  %                       |          |          |
  %                       |----------|----------|
  %                       |          |          |
  %                       |   p=3    |   p=4    |
  %                       |          |          |
  %                       |----------|----------|
  %
  % You can use sublot(2,2,1) or subplot(221) to activate the frames.
  
  figure(1)
  % Plot in the first window p=1
  subplot(2,2,1);
  imagesc(echo_data.GPS_time,echo_data.Depth,db(echo_data.Data,'power'))
  title('Echogram');
  % Plot in the second window p=2
  subplot(2,2,2);
  plot(echo_data.GPS_time)
  title('GPS Time')
  % Plot in the third window p=3
  subplot(2,2,3);
  plot(echo_data.GPS_time,echo_data.Surface,'k-')
  title('Surface')
  % Plot in the fourth window p=4
  subplot(2,2,4);
  plot(echo_data.GPS_time,echo_data.Bottom,'k-')
  title('Bottom')
  
  close all;
  
  
  % PLOTTING CRESIS LAYER DATA
  % =======================================================================
  
  % Say we want to plot the surface, and bottom on a single
  % figure with a common axes.
  
  file_path = fullfile(data_folder_path,'CRESIS_CSV_DATA.csv');
  fid = fopen(file_path);
  csv_data = textscan(fid,'%f%f%f%f%f%f%f%f%s','delimiter',',','headerlines',1);
  fclose(fid);
  
  % Note the colums of a CReSIS CSV File
  % LAT,LON,TIME,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY
  lat = csv_data{1};
  lon = csv_data{2};
  utc_time = csv_data{3};
  thickness = csv_data{4};
  elev = csv_data{5};
  frame = csv_data{6};
  surface = csv_data{7};
  bottom = csv_data{8};
  quality = csv_data{9};
  % This is a faster way to do the same thing:
  [lat,lon,utc_time,thickness,elev,frame,surface,bottom,quality] = deal(csv_data);
  % Another way to do it:
  [lat,lon,utc_time,thickness,elev,frame,surface,bottom,quality] = csv_data{:};
  
  % Now we can plot surface and bottom. It's important to note in CReSIS data
  % variables surface and bottom are the distance from elevation to each
  % respective surface. We need to calculate the actual values.  
  a_surf = csv_data{5}-csv_data{7};
  a_bed = csv_data{5}-csv_data{8};

  figure(1);
  plot(csv_data{3},a_surf,'r-')
  hold on;
  plot(csv_data{3},a_bed,'g-')
  
  % Lets put a legend,title,and labels on the plot.
  title('Surface and Bottom Layers');
  xlabel('UTC Time (SOD)');
  ylabel('Elevation WGS84 (m)');
  legend('Surface','Bottom');   % Go in the order of the plots.
  grid on;
  xlim(csv_data{3}([1 end])); % Set the x-limits to fit the data
  ylim([min(min(a_surf,a_bed))-100 max(max(a_surf,a_bed))+100]); % Set the y-limits with the 100 m buffer
  
  close all;
  
  
  % STATISTICAL PLOTS
  % =======================================================================
  load seamount; % This file contains topographical data for a seamount
  fprintf('\n')  % and it ships with the MATLAB package.
  scatter(x,y,5,z)              % Creates a scatter plot. The third dimension is
  % represented by color. Since the z vector
  % has the same dimensions as the x and y
  % vectors, color is mapped linearly over the
  % scatter plot.
  figure(2);
  x = -2.9:0.2:2.9;    % Makes a bar graph showing a two-sided exponential function.
  bar(x,exp(-x.*x),'r')% The 'r' parameter sets its color to red and the x vector
  % specifies the spacing. We determine the Y
  % values inside of the bar() call rather than
  % before the call as we have done previously.
  
  % Plot histogram of discrete uniform distribution from 1 to 6
  figure(3);
  x = 1:6;
  y = ceil(rand(1,10000)*6);
  hist(y,x);
  
  close all;

  
  % SIMPLE INTERACTIVE PLOT
  % =======================================================================
  
  % The ginput() function allows a user to draw input to a figure which
  % matlab will store in a variable.
  
  % The following is a simple example that will allow a user to draw 5
  % points on a new figure.
  
  figure(1);                % Create a new figure
  fprintf('Input 5 points by clicking on the figure.\n');
  [X Y] = ginput(5);        % Get input 5 times, store in X,Y
  plot(X,Y,'r+')            % Plot the input with a red +
  
  help impoly
  % Another common input tool is the polygon tool. Read the help for
  % details. Left click to places vertices. Double click to add the last
  % vertex.
  H = impoly
  
  close all;
  
  
  % SAVE FIGURES
  % ---------------------------------------------------------
  
  % See the following for specifics on saving images in MATLAB
  % doc print
  % doc saveas
  
  C = load('clown');    % Loads the MATLAB file
  Cl = image(C.X);
  
  % Save as JPEG
  out_fn = fullfile(data_folder_path,'clown_image')
  saveas(Cl,out_fn,'jpeg');
  
  % Save as Matlab .fig file (contains all the information to recreate, but
  % requires Matlab to load)
  out_fn = fullfile(data_folder_path,'clown_image')
  saveas(Cl,out_fn,'fig');
  
  close all;
  
  
  % 3D PLOTS
  % ---------------------------------------------------------
  % See doc plot3 for more information.
  
  fprintf('Plot a 3D Helix.\n');
  
  figure(1);              % Plot a helix. With 3D plots, the MATLAB figure
  t = 0:pi/50:10*pi;      % window tools are very useful: use the Hand tool
  plot3(sin(t),cos(t),t)  % to move the plot around, the Rotate tool to view it
  axis square;            % from different angles, and the Zoom tool to adjust
  grid on                 % the plot scale. They are located right below the menu
  % toolbar.
  
  
  % SURFACE PLOTS
  % ---------------------------------------------------------
  % See doc mesh; doc surf; doc contour;
  
  fprintf('Plot a 3D Mesh and 2D contours.\n');
  
  figure(1);
  [X,Y] = meshgrid(-3:.125:3);    % Create evenly-spaced X and Y grids
  Z = peaks(X,Y);                 % Generate a Z grid using the peaks function
  meshc(X,Y,Z)                    % Create a mesh plot with a contour map below
  axis([-3 3 -3 3 -10 5])
  
  fprintf('Plot a 3D Mesh and 2D curtain.\n');
  
  [X,Y] = meshgrid(-3:.125:3);
  Z = peaks(X,Y);
  meshz(X,Y,Z)                    % Same as above, but now a curtain is added instead of a contour map
  
  colormap hsv                    % Set hsv colormap
  axis([-3 3 -3 3 -10 5])
  
  fprintf('Plot a 2D peaks contour with labels.\n');
  
  Z = peaks;                      % Get default size (49x49) Z grid from peaks function
  [C,h] = contour(interp2(Z,4));  % Generate a contour plot for these Z values... the interp2
  % function is a 2D interpolation that
  % basically expands the number of values
  % in our Z grid by guessing what value(s)
  % would naturally occur between adjacent
  % matrix elements. So you can think of it
  % as increasing the precision of our peaks
  % function so that our contour map is
  % more precise
  
  clabel(C,h);      % Turn on contour labels so that the levels
  % of each contour are described
  % quantitatively
  
  fprintf('Section 5 Complete.\n');
  return
end

%% Section 6: Mapping toolbox
if section_number == 6
  close all
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 6: Mapping toolbox\n');
  fprintf('========================================================\n')
  
  keyboard

  % =======================================================================
  fprintf('Load and Plot World Coastline Data.\n');
  
  load coast;         % Loads built in coastline vector data
  
  axesm mercator;     % Sets the projection to mercator
  framem;             % Adds a frame to the map window
  plotm(lat,long);    % Plots the lat/long of the coastline data.
  
  
  fprintf(' Clean Up:\n')
  clear lat long; close all;
  
  
  
  fprintf('Load the Mapping GUI (See code comments).\n');
  
  load coast;         % Loads build in coastline vector data
  axesm;              % axesm by itself loads the Mapping GUI
  
  % Explore the dialog options of the Mapping GUI
  
  % ===== TEST =====
  % Create a map, using the mapping GUI that contains the following features:
  % 1) Uses the Mercator Cylindrical projection
  % 2) Has a simple frame
  % 3) Has a grid with dashed lines.
  %     - 60 deg meridians and 30 deg parallels.
  % 4) Has labels.
  %     - 60 deg meridians and 30 deg parallels.
  % 5) Shows the coastline vector data
  
  % If you cannot get your map to work correctly, or you are not sure what it
  % should look like use the code below to generate the map, then close and
  % re-try creation with the GUI.
  
  
  % =======================================================================
  fprintf('Load the correct map w/o the GUI.\n');
  
  load coast;
  axesm('MapProjection','mercator','Frame','on','Grid','on',...
    'MLineLocation',60,'PLineLocation',30,'MeridianLabel','on',...
    'ParallelLabel','on','MLabelLocation',60,'PLabelLocation',30);
  plotm(lat,long);
  
  % Clean Up:
  clear lat long; close all;
  
  
  % =======================================================================
  fprintf('Load the korea grid and plot a texture map.\n');
  
  load korea;             % Load the korea topographic grid
  figure;                 % Create a new figure window
  worldmap(map,refvec);   % Create a worldmap reference map.
  
  % Note "refvec" is part of the file "korea". It is the reference vector
  % used to display the topographic and bathymetric data.
  
  % Display the data as a "Texture Map"
  geoshow(gca,map,refvec,'DisplayType','texturemap');
  
  % Note: "gca" is a pointer to the current figure window. It basically tells
  % the data where to go.
  
  colormap(demcmap(map)); % Sets a terrain appropriate colormap.
  
 
  % Clean Up:
  clear description map maplegend refvec source; close all;
  
  
  % =======================================================================
  fprintf('Load and plot a map with composite data types.\n');
  
  % First we will re-load and plot the Korea Map.
  
  load korea;             % Load the korea topographic grid
  figure;                 % Create a new figure window
  worldmap(map,refvec);   % Create a worldmap reference map.
  % Display the data as a "Texture Map"
  geoshow(gca,map,refvec,'DisplayType','texturemap');
  colormap(demcmap(map)); % Sets a terrain appropriate colormap.
  
  % Now, Lets read a shapefile called "landareas" and plot the boundaries as
  % black lines on top of the korea plot.
  
  landareas = shaperead('landareas','UseGeoCoords',true);   % Load the shapefile.
  
  % First lets explore the shapefile we just loaded. The file is in an array
  % structure. Lets access the first feature and see whats included.
  
  landareas(1)  % Display the properties of the first feature.
  
  % Note the feature has a type (Geometry), BoundingBox, Lon & Lat vector and
  % a Name attribute. This is the basic "Geographic Data Structure"
  
  % Plot the landareas in black using geoshow.
  geoshow([landareas.Lat],[landareas.Lon],'Color','black');
  
  
  % Clean Up:
  clear landareas description map maplegend refvec source ans; close all;
  
  
  % =======================================================================
  fprintf('Load and plot a 3D map.\n');
  
  % The DEM data is provided with MATLAB in ZIP file format. We first need to
  % unzip and store the data.
  
  % Unzip the file, and store the files it contains in "filenames"
  filenames = gunzip('sanfranciscos.dem.gz', tempdir);
  demFilename = filenames{1};   % Select the filename we want, store in demFilename.
  
  % The data is a USGS 1:24,000 DEM File, MATLAB has a built in tool called
  % usgs24kdem.m that reads this format. Lets read the file.
  
  [lat,lon,elev] = usgs24kdem(demFilename,1);
  
  % Note the "2" at the end of the function call. This is a sampling factor.
  % Since it is set to 2 we will actually ready every other point from the
  % file. This will save us some time.
  
  % Lets delete the temporary file we unzipped.
  delete(demFilename);
  
  % To help visualize in 3D lets set all elevation values at sea level to -1
  % (This will become "blue" in the colormap.)
  elev(elev==0) = -1;
  
  % We next need to compute the latitude and longitude limits of the DEM.
  latlim = [min(lat(:)) max(lat(:))];
  lonlim = [min(lon(:)) max(lon(:))];
  
  % We can now begin the 3D plot.
  figure;        % Open a new figure window.
  
  % Function usamap.m creates an empty map axes (same as axesm) with specific
  % projection and location information. See "help usamap" for more info.
  usamap(latlim,lonlim);
  
  geoshow(lat,lon,elev,'DisplayType','Surface');  % Plot the elvation at lat,lon
  demcmap(elev); % Set the elevation colormap
  
  % Next we need to set the vertical exageration of the DEM.
  daspectm('m',2); % Exageration in meters by a factor of 2
  view(3);    % Finally we can set the view.
  
  % Lets add a colorbar to the map
  colorbar;   % It does what it is, or is what is does, or is it?
  

  % Clean Up:
  clear demFilename elev filenames lat latlim lon lonlim; close all;
  
  
  % =======================================================================
  fprintf('Read and explore and ESRI Shapefile.\n');
  
  % There are two main shapefile loading functions.
  % (1) shapeinfo.m
  % (2) shaperead.m
  %
  % shapeinfo loads the information , not the data, of a shapefile.
  
  % Load the info of the "concord_roads" dataset and explore the structure.
  
  info = shapeinfo('concord_roads');
  
  % Info is a structure with the following fields:
  % (1) Filename: Character array of the filenames associated with the
  %               shapefile.
  % (2) ShapeType: Type of feature (Point,Polygon,PolyLine,etc...)
  % (3) BoundingBox: Array of the bounding box coordinates.
  %                  [ min(X) min(Y) ]
  %                  [ max(X) max(Y) ]
  % (4) NumFeatures: The number of features in the shapefile.
  % (5) Attributes: Structure Array with the name and type of each attribute.
  
  % Explore the information and fields stored in the info structure. Then
  % continue.
  
  % Load the data of the "concord_roads" dataset and explore the structure.
  
  roads = shaperead('concord_roads');   % Loads the concord_roads shapefile.
  close all;
  mapshow(roads);                       % Since we have XY coordinates (Not Lat/Lon) we use mapshow)
  
  % Roads is a structure array with a structure for every feature.
  % Note that length(roads) == info.NumFeatures
  
  % Lets take a look at a single feature of the roads dataset.
  road10 = roads(10);
  
  % Road 10 Is a segment of the "Wright Farm" road. A single structure of the
  % shapefile contains the following fields:
  
  % DEFAULT (WITH EVERY DATASET)
  % (1) Geometry: Type of feature (Point,Polygon,PolyLine,etc...)
  % (2) BoundingBox: Array of the bounding box coordinates.
  %                  [ min(X) min(Y) ]
  %                  [ max(X) max(Y) ]
  % (3) X: X Coordinate/s of Feature
  % (4) Y: Y Coordinate/s of Feature
  
  % OPTIONAL (VARY BY DATASET)
  % There will be a single field for every additional attribute in the
  % dataset. Ex. (STREETNAME, RT_NUMBER, CLASS, ADMIN_TYPE, LENGTH, etc...)
  
  % What if the shapefile we are loading contains Lat\Lon attributes and not
  % X\Y coordinates. How do we use those?
  
  world_cities = shaperead('worldcities', 'UseGeoCoords', true);
  
  % Lets plot the world_cities on a world land areas map using geoshow.
  geoshow('landareas.shp');   % This will load and plot the land areas shapefile.
  geoshow(world_cities);    % This will plot the cities shapefile. Note since we have Lat\Lon we use geoshow.m
  

  % Clean Up
  clear info roads10 roads world_cities; close all;
  
  
  % =======================================================================
  fprintf('Read and explore a GeoTIFF.\n');
  
  % There are two main GeoTIFF loading functions.
  % (1) geotiffinfo.m
  % (2) geotiffread.m
  
  % Lets get the info for the "boston.tif" geotiff.
  
  info = geotiffinfo('boston.tif');
  
  % Take a look at info, you will notice quite a few fields with lots of
  % information. This information will vary by dataset.
  
  % Lets load the "boston.tif" geotiff.
  [boston R bbox] = geotiffread('boston.tif'); % Loads the data into "boston" the reference matrix into "R" and the bounding box into "bbox"
  
  % Lets plot the boston geotiff with mapshow.
  figure;               % Create a new figure;
  mapshow(boston,R);    % Plot the data
  axis image off        % Turn the axis to image mode and then off

  clear R bbox boston info; close all;
  
  
  % =======================================================================
  fprintf('Load a GeoTIFF and Shapefile, plot the composite.\n');
  
  % Lets load the boston roads shapefile.
  boston_roads = shaperead('boston_roads');
  
  % Lets plot the geotiff and then roads on top.
  figure;
  mapshow boston.tif;
  axis image off;
  
  % You should from geotiffinfo and shapeinfo that the coordinates are
  % different for the files. The tiff is in feet and the shapefile in meters.
  % We need to convert the roads to overlay the two.
  
  Ft_to_M = unitsratio('sf','meter');   % Get the conversion factor from foot to meter
  x = Ft_to_M * [boston_roads.X];       % Convert roads X (in feet) to x (in meters)
  y = Ft_to_M * [boston_roads.Y];       % Convert roads Y (in feet) to y (in meters)
  
  mapshow(x,y);   % Plot the new roads over the GeoTIFF.

  % Clean Up
  clear Ft_to_M boston_roads x y; close all;
  
  
  % =======================================================================
  fprintf('Lets write some XY data to an ESRI shapefile...\n');
  
  % Run the following commands to set up a dataset to save.
  st = shaperead('usastatehi');
  xk = [st(16).X]; xk = xk(~isnan(xk)); xk = xk(1:end-1);
  yk = [st(16).Y]; yk = yk(~isnan(yk)); yk = yk(1:end-1);
  xm = [st(25).X]; xm = xm(~isnan(xm)); xm = xm(1:end-1);
  ym = [st(25).Y]; ym = ym(~isnan(ym)); ym = ym(1:end-1);
  clear st;
  
  % Great! So were starting with variables xk,yk and xm,ym. Strings of coordinates
  % representing the polygon boundary of Kansas and Missouri.
  
  % Plot xc and yc to confirm the data. You should see Kansas and Missouri
  % with a gap in the NW Corner of Kansas.
  plot(xk,yk); hold on; plot(xm,ym)
  
  % To convert from vectors of XY coordinates to a SHAPEFILE a few conditions
  % must be met.
  % (1) Starting coordinate and ending coordinate must be the same.
  % (2) The end value of a polygon should be a NaN
  
  % We first should check if the starting values equal the ending values.
  
  xk(1) == xk(end);   % If this equals 1 the coordinates match, if 0 they dont.
  xm(1) == xm(end);   % If this equals 1 the coordinates match, if 0 they dont.
  
  % Notice neither meet the first condition. Lets fix this. We need to add a
  % new value (end+1) and set it to the first value (1).
  xk(end+1) = xk(1);
  yk(end+1) = yk(1);
  xm(end+1) = xm(1);
  ym(end+1) = ym(1);
  
  % We next need to add the NaN to the end of each coordinate string. Having
  % NaN at the end of a coordinate string represents the end of a polygon.
  xk(end+1) = NaN;
  yk(end+1) = NaN;
  xm(end+1) = NaN;
  ym(end+1) = NaN;
  
  % Next we need to construct the Geographic Data Structure. Each polygon
  % must contain the following fields within the structure:
  % (1) "Geometry": One of the following shape types: 'Point', 'MultiPoint', 'Line', or 'Polygon'.
  % (2) "BoundingBox": Specifies the minimum and maximum feature coordinate
  %                    values in each dimension in the following form: [min(X) min(Y) ; max(X) max(Y)]
  % (3) "X/Y" or "LON/LAT": 1-by-N coordinate vector.
  % (4) "Attributes": Any additional attribute data would go here.
  
  % Note a mapstruct and geostruct are the same. A mapstruct contains XY and
  % a geostruct contains Lon/Lat.
  
  % Lets construct a 2x1 array with fields
  % Geometry,BoundingBox,X,Y,STATENAME.
  
  mapstruct(1,1).Geometry = 'Polygon'; % Set both Geometry's to Polygon
  mapstruct(2,1).Geometry = 'Polygon';
  
  % Lets find the bounding coordinates for Kansas and Missouri.
  k_bbox = [min(xk) min(yk) ; max(xk) max(yk)];
  m_bbox = [min(xm) min(ym) ; max(xm) max(ym)];
  
  % Now we can add BoundingBox to the 2x1 array.
  mapstruct(1,1).BoundingBox = k_bbox; % We'll make Kansas #1.
  mapstruct(2,1).BoundingBox = m_bbox;
  
  % Next lets add the actual coordinate strings to the 2x1 array.
  mapstruct(1,1).X = xk;  mapstruct(1,1).Y = yk;
  mapstruct(2,1).X = xm;  mapstruct(2,1).Y = ym;
  
  % Finally we can add the Name Attribute.
  mapstruct(1,1).Name = 'Kansas';
  mapstruct(2,1).Name = 'Missouri';
  
  % At this point the map structure should be ready for writing. If you type
  % "mapstruct" in the command line you should get:
  %   mapstruct =
  %   2x1 struct array with fields:
  %       Geometry
  %       BoundingBox
  %       X
  %       Y
  %       Name
  
  % If thats correct, lets write the file!  
  % First specify and output location (Change the path) and filename:
  out_shp = fullfile(data_folder_path,'out_NewSHP');
  
  % We can now use the shapewrite.m function to write the new shapefile.
  shapewrite(mapstruct,out_shp);
  
  % Now lets load the shapefile we just wrote and plot the data using Mapshow.
  ksmo_bound = shaperead(out_shp);
  mapshow(ksmo_bound,'DisplayType','polygon');
  
  % Examine the ksmo_bound file to see how the mapstruct we created is
  % converted in a shapefile, then clean up! Make sure to delete the
  % shapefile from whatever directory you created it in "out_shp"
  
  % Clean Up
  clear ans ax k_bbox ksmo_bound m_bbox mapstruct out_shp xk xm yk ym; close all;
  
  
  % =======================================================================
  fprintf('Lets do some coordinate conversion ...\n');
  
  % CReSIS has a few conversion functions in the CReSIS-Toolbox.
  % These functions are:
  % (1) geodetic_to_utm.m
  % (2) geodetic_to_stereographic.m
  % (3) geodetic_to_along_track.m
  
  % These functions are used to convert from geodetic (Geographic)
  % coordinates in the WGS1984 datum to projected coordinates or along-track
  % distances. We'll go through each function individually.
  
  % Let's load the Antarctica Coastline and convert the coordinates to
  % sterographic. Specifically UPS (Universal Polar Stereographic) South.
  
  coastline = shaperead('landareas', 'UseGeoCoords', true,...
    'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});
  
  % If we inspect the coasline data both visually and by looking at the data
  % we see that the Lat/Lon coordinates are Geographic (Geodetic).
  
  geoshow(coastline);
  
  inspect_lat = coastline.Lat(1);
  inspect_lon = coastline.Lon(1);
  
  % This is not a good way to visualize or work with Antarctica. We should
  % convert the coordinates to UPS South coordinates.
  
  [x,y,mstruct] = geodetic_to_stereographic(coastline.Lat,coastline.Lon);
  
  % To visualize the change use this simple plot function.
  figure;
  subplot(2,1,1); plot(coastline.Lon,coastline.Lat); title('Geodetic Coordinates');
  subplot(2,1,2); plot(x,y); title('Stereographic Coordinates');
  
  % How does this work? Below is the entire geodetic_to_stereographic.m
  
  % ===================================
  % mstruct = defaultm('ups');
  % if mean(lat(:)) > 0
  %   mstruct.zone = 'north';
  % else
  %   mstruct.zone = 'south';
  % end
  % mstruct = defaultm(mstruct);
  % mstruct.falsenorthing = 0;
  % mstruct.falseeasting = 0;
  % physical_constants;
  % mstruct.geoid = WGS84.ellipsoid;
  %
  % [x,y] = mfwdtran(mstruct,lat,lon);
  % ===================================
  
  % As you can see there is actual only one function call to mfwdtran.m and
  % the rest is just the setup of a mapstruct. MFWDTRAN simple takes the
  % mapstruct projection information and transforms the latitude and longitude
  % into that coordinate system. So geodetic_to_stereographic.m and
  % geodetic_to_utm.m are just specific functions designed using MFWDTRAN.
  
  % geodetic_to_along_track.m is a a bit different in the fact that it does
  % not convert to projected coordinates rather to a vector of linear
  % distances between coordinates "along a track".
  
  
  % =======================================================================
  fprintf('Lets use the geodetic_to_along_track function ...\n');
  
  % Say we have the following two coordinate strings and want to find the
  % distance between each pair of coordinates.
  
  lat = [-83.8906  -84.4575  -84.4119  -84.6315 -84.7173  -84.8966  -85.0214];
  lon = [-179.9003 -178.2340 -175.0965 -169.0558 -167.5598 -166.9709 -163.5254];
  elev_ft = [439 678 260 273 273 265 266];
  
  % First lets plot to visualize this problem.
  figure;
  plot(lon,lat,'b-'); hold on; plot(lon,lat,'k+');
  
  % So we want the length of the blue line between each +. To find this we
  % use the geodetic_to_along_track.m function. The function needs lat,lon
  % and elevation in meters. So we need to convert out elevation.
  
  elev_m = elev_ft*unitsratio('ft','m');
  
  along_track_dist = geodetic_to_along_track(lat,lon,elev_m);
  
  % So if we look at along_track_dist now we see that the result is the
  % cummulative distance, not the difference between each point.
  
  % along_track_dist =
  %    1.0e+05 *
  %          0    0.6609    1.0048    1.6938    1.8761    2.0850    2.4507
  
  % Since we want the distance between each point we can use the diff.m
  % function to get the difference between each vale.
  
  dist_between_coords = diff(along_track_dist);
  
  % dist_between_coords =
  %    1.0e+04 *
  %     6.6088    3.4395    6.8892    1.8232    2.0894    3.6573
  
  % With that complete you should have a fundemental understanding of
  % projecting coordinates and converting coordinates to cummulative and
  % between point distances using CReSIS functions.
  
  % Clean Up
  clear along_track_dist coastline dist_between_coords elev_ft elev_m inspect_lat inspect_lon lat lon mstruct x y; close all;
  
  
  % =======================================================================
  fprintf('Lets do a coordinate projection ...\n');
  
  % MATLAB has a few important built in coordinate conversion and projection
  % functions that are important to note. Thos are:
  % (1) projfwd.m
  % (2) projinv.m
  % (3) mfwdtran.m
  % (4) minvtran.m
  
  % These 4 functions break down into two groups:
  % (1) Coordinate Projection
  % (2) Coordinate Transformation
  
  % Coordinate Projection (projfwd,projinv)
  % -------------------------------------------------------------------------
  % Say we have the following WGS84 Coordinate pair.
  
  lat = 38.9521;
  lon = -95.2640;
  
  % To view a complete list of the possible projections:
  projlist;
  
  % Project into the Transverse Mercator projection.
  proj = defaultm('tranmerc');
  [x,y] = projfwd(proj,lat,lon);
  
  % Great now we have the follwoing coordinate pair:
  % x = -1.0312 , y = 1.6838
  
  % What is we started with the Transverse Mercator coordinates and wanted
  % LAT/LON in degrees, WGS1984?
  
  [lat2,lon2] = projinv(proj,x,y);
  
  % So projinv is simply the inverse of projfwd. You can confirm that lat2 =
  % lat and lon2 = lon if you dont believe it!
  
  % Clean Up
  clear lat lon x y lat2 lon2 proj;
  
  
  fprintf('Lets do a coordinate transformation ...\n');
  
  
  % Coordinate Transformation (mfwdtran,minvtran)
  % =======================================================================
  % Say we have the following WGS84 Coordinate pair.
  
  lat = 38.9521;
  lon = -95.2640;
  
  % Transform these cordinates into UTM map coordinates. First we need
  % to set up our map projection structure.
  
  mstruct = defaultm('utm');          % Load the default UTM mstruct.
  mstruct.zone = utmzone(lat,lon);    % Set the UTM Zone.
  mstruct = defaultm(mstruct);        % Re-Set the mstruct with the new zone.
  
  [x,y] = mfwdtran(mstruct,lat,lon);
  
  % Great! We now have the following coordinate pair:
  % x = 303810 , y = 4313700 , zone 15S
  
  % Just like before what if we want to convert from x,y back to lat,lon?
  
  [lat2,lon2] = minvtran(mstruct,x,y);
  
  % Just like before minvtran is the inverse of mfwdtran. Compare lat2/lat
  % and lon2/lon to confirm!
  
  % Clean Up
  clear lat lon x y lat2 lon2 mstruct;
  
  
  fprintf('Lets apply all of the above to CReSIS data ...\n');
  
  if ispc
    csv_path = 'X:/public/data/rds/2011_Antarctica_DC8/csv/Browse_2011_Antarctica_DC8.csv';
  else
    csv_path = '/cresis/snfs1/public/data/rds/2011_Antarctica_DC8/csv/Browse_2011_Antarctica_DC8.csv';
  end
  gtiff_path = ct_filename_gis([],fullfile('antarctica','Landsat-7','Antarctica_LIMA_peninsula.tif'));
  
  if exist(csv_path,'file') && exist(gtiff_path,'file')
    % Assign a file identifier and open a CSV file of flight lines as read-only
    fid = fopen(csv_path,'r');
    
    % Create variable csv_file that reads the csv in based upon data type (%f
    % is floating point and %s is a string), that the file has a headerline at
    % the top, and that the values are delimited (i.e. separated) by commas.
    csv_file = textscan(fid,'%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
    
    lat = (csv_file{1}(1:5:end)); % Assign all values from column 1 as "lat" & only take every 5th point (speeds up load time)
    lon = (csv_file{2}(1:5:end)); % Assign all values from column 2 as "lon" & only take every 5th point (speeds up load time)
    
    fclose(fid); % Close CSV file
    
    % Grab the GeoTIFF's info and store in in structure proj
    proj = geotiffinfo(gtiff_path);
    
    figure(1);  % All plots will appear in this figure window
    hold on;    % Forces all other plot commands to stay on same figure window
    RGB = {};
    R = {};
    [RGB, R, tmp] = geotiffread(proj.Filename); % Reads the projection information into a variable that can be plotted
    R = R/1e3;  % Convert projection data from meters to kilometers
    
    mapshow(RGB,R); % Plots the GeoTIFF
    
    % Adjust the axis so it reads in kilometers
    axis([R(3,1)+[0 R(2,1)*(size(RGB,2)-1)] sort(R(3,2)+[0 R(1,2)*(size(RGB,1)-1)])]);
    
    % Apply projection data 'proj' and name the new variable set x and y
    [x,y] = projfwd(proj,lat,lon);
    
    % Adjust the x/y values so they can be plotted on a grid set for
    % kilometers
    x = x/1e3;
    y = y/1e3;
    
    flight_lines = plot(x,y,'r-'); % Plot projected coordinate pair on to the figure
    
    % Add a title, x and y-axis labels to the figure
    map_title = title('2009 Antarctica TO Browse File');
    
    % Add axis labels
    map_x = xlabel('X (km)');
    map_y = ylabel('Y (km)');
    
    % Add a legend for the flight lines, and place it in the most ideal
    % location ('Location','Best')
    flight_legend = legend([flight_lines],'Flight Lines','Location','Best');
    
    % Change figure position and size to show up correctly in plot window
    fig_position = [0.5 2.5 6 8];
    fig_size = [50,50,600,800];
    set(figure(1),'PaperPosition',fig_position);
    set(figure(1),'Position',fig_size);
    
    % Output figure(1) as a jpeg with a resolution of 450 DPI
    output_path_and_file = fullfile(data_folder_path,'newfile.jpg');
    print(figure(1),'-djpeg','-r450',output_path_and_file);
  end
  
  fprintf('Section 6 Complete.\n');
  return
end

%% Section 7: Vectorization
if section_number == 7
  close all
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 7: Vectorization \n');
  fprintf('========================================================\n')
  
  fprintf(' See the following video and page for detailed information about improving\n')
  fprintf(' performance in MATLAB. This is important as you begin to work with\n')
  fprintf(' large datasets.\n')
  fprintf('\n');
  fprintf('http://www.mathworks.com/company/newsletters/articles/programming-patterns-maximizing-code-performance-by-optimizing-memory-access.html\n')
  fprintf('\n');
  fprintf('This section will explain some methods of optimization using MATLAB.\n')
  
  keyboard
  
  % PREALLOCATION AND VECTORIZATION 1
  % =======================================================================
  % This example shows the value of preallocating a matrix before filling
  % it in one element at a time. It also shows the importance of
  % vectorization.
  %
  % Each of the three sets of operations creates an identical sine wave.
  % Look at how this is done and compare the computation time for each.
  
  fprintf('\nRun "dbcont" to run the following for loops.\n');
  
  fprintf('\nFor loop with no preallocation:\n');
  
  clear x1 x2 x3;
  len  = 1e6;
  tic
  for idx0 = 1:len
    x1(idx0) = cos(2*pi*100*idx0); % Data is created and stored on the fly.
  end
  toc
  
  fprintf('For loop with preallocation:\n');
  
  tic
  x2 = zeros(1,len);    % We have pre-allocated a set of zeros to hold the data.
  for idx0 = 1:len
    x2(idx0) = cos(2*pi*100*idx0);
  end
  toc
  
  fprintf('Vectorized function:\n');
  
  tic
  t = 1:len;
  x3 = cos(2*pi*100*t);   % Here we use the cos() function. This function is vectorized.
  toc
  
  keyboard
  
  
  % PREALLOCATION AND VECTORIZATION 2
  % =======================================================================
  % This example shows the value of preallocating a matrix before filling
  % it in one element at a time. It also shows the importance of
  % vectorization.
  %
  % Each of the three sets of operations computes an element wise matrix
  % multiply. Look at how this is done and compare the computation time for
  % each.
  
  fprintf('\nRun "dbcont" to run the following for loops.\n');
  
  len = 100;
  
  clear x1 x2 x3;
  r1 = randn(len,len); % Create len x len matrix of Gaussian random variables
  r2 = randn(len,len); % Create len x len matrix of Gaussian random variables
  
  fprintf('\nFor loop without preallocation.\n')
  tic
  for idx0 = 1:len
    for idx1 = 1:len
      x1(idx0,idx1) = r1(idx0,idx1)*r2(idx0,idx1);
    end
  end
  toc
  
  fprintf('For loop with preallocation.\n')
  tic
  x2 = zeros(len);
  for idx0 = 1:len
    for idx1 = 1:len
      x2(idx0,idx1) = r1(idx0,idx1)*r2(idx0,idx1);
    end
  end
  toc
  
  fprintf('Vectorized.\n')
  tic
  x3 = r1.*r2;
  toc
  
  keyboard
  
  
  % PERMUTE
  % =======================================================================
  
  % The physical memory on the computer is best thought of as a single long
  % vector starting at address zero and going up to an address that is as
  % big as the memory you have.  Matrices are stored in memory by mapping
  % each dimension out into one long vector. For example:
  A = [1 2 3; 4 5 6; 7 8 9]
  % is stored in memory like this:
  A(:)
  % Rows are traversed first, then columns, then each upper dimension in
  % turn.
  
  % Indexing a matrix can either use the more intuitive (row,column) form
  % or the vectorized form. The vectorized form is used any time the
  % dimensions of the indexing are less than the dimensions of the matrix.
  % Here are two ways to access the same element:
  A(1,2)
  A(4)
  % Can you determine a general code for converting between the two
  % indexing methods?
  
  % Computers use a tiered memory system:
  %   L1 Cache: very fast, no delay to use variables here
  %   L2 Cache: fast, small delay to use variables here
  %   RAM: fast, big delay to use variables here
  %   Hard Drive: slow, big delay to use variables here
  % The fast the tier, the less memory that is available. Computers use
  % algorithms based on locality in space and time (i.e. variables that
  % are close in memory are more likely to be used close in time) to decide
  % which variables should be placed in the faster caches at any given
  % time.
  %
  % Since matrices are vectorized along rows first (i.e. row-major), when
  % a single element of a row is accessed it will likely grab the entire
  % row into memory. This is really helpful if you are doing operations
  % along rows, but not helpful if you are doing operations along columns.
  
  fprintf('\nRun "dbcont" to run the following operations.\n');
  
  % Row-wise FFT is faster
  clear A B;
  A = randn(1e4,1e4);
  tic;
  B = fft(A,[],1); % FFT along first-dimension
  toc;
  
  % Column-wise FFT is slower
  clear A B;
  A = randn(1e4,1e4);
  tic;
  B = fft(A,[],2); % FFT along second-dimension
  toc;
  
  keyboard
  % BSXFUN (Binary Singleton Function)
  % =======================================================================
  % This is commonly used to apply element wise operations on every
  % possibly combination of the elements from two vectors.
  
  % Syntax: A=bsxfun(fun, x, y). where A, x, and y are vectors or matrixes.
  R = [10 11 12].'
  T = [1 2 3]
  C = bsxfun(@minus, R, T) % Performs R minus T for every combination of R and T
  
  % This is the same as replicating the vectors first to match in size
  % and then applying the operation:
  repmat(R,[1 3])
  repmat(T,[3 1])
  C_equivalent = repmat(R,[1 3]) - repmat(T,[3 1])
  
  % Another way of using Bsxfun is to specify an "anonymous funtion". This
  % lets you create your own custom function.
  R = [0 1 2].';
  T = [1 2 3];
  C = bsxfun(@(R,T) 2*R+T, R,T)
  
  repmat(R,[1 3])
  repmat(T,[3 1])
  C_equivalent = 2*repmat(R,[1 3]) + repmat(T,[3 1])
  
  % Example which calculates the surface clutter angle for altitude above
  % groud level H and ice thickness T assuming a flat surface.
  Height = [200:1000];
  Thickness = [500:2000];
  physical_constants;
  speed_in_ice=c/sqrt(er_ice);
  ca = bsxfun(@(Height,Thickness) acos(Height ./ ((Height/c+Thickness/speed_in_ice)*c)), Height, Thickness.');
  imagesc(Height, Thickness, ca*180/pi)
  xlabel('Height (m)');
  ylabel('Thickness (m)');
  hcolor = colorbar;
  set( get(hcolor,'YLabel') , 'String', 'Clutter angle (deg)');
  
  keyboard
  
  % REPMAT, BSXFUN to vectorize commands
  % =======================================================================
  % Commonly one variable of an operation is a vector and another variable
  % is a matrix. For example: Let A be a matrix that we want to remove
  % the mean of each row from.
  A = randn(10000,1000);
  row_mean = mean(A,2);
  size(A)
  size(row_mean)
  
  % We cannot do A - row_mean because the dimensions to not match. There
  % are a couple of common ways to do this in a vectorized way:
  fprintf('\nRun "dbcont" to run the following operations.\n');
  
  clear A1 A2 A3;
  % 1. Use replicate matrix to replicate row_mean
  tic
  A1 = A - repmat(row_mean,[1 size(A,2)]);
  toc
  % 2. Use bsxfun
  tic
  A2 = bsxfun(@minus,A,row_mean);
  toc
  % 3. For loop
  tic;
  A3 = zeros(size(A)); % Preallocate
  for row = 1:size(A,1)
    % Subtract the mean one row at a time (this works because now we are
    % subtracting a scalar)
    A3(row,:) = A(row,:) - row_mean(row);
  end
  toc;
  
  fprintf('Section 7 Complete.\n');
  return
end

%% Section 8: Debugging
if section_number == 8
  
   fprintf('\n=========================================================\n');
  fprintf(' Section 8: Debugging\n');
  fprintf('========================================================\n')
  
  fprintf(' Debugging tutorial:\n')
  fprintf(' http://www.youtube.com/watch?v=Z4vFymKhNno\n')
  fprintf('\n')
  fprintf(' There are a few other useful things to note here:\n')
  fprintf('\n')
  fprintf(' TRY/CATCH Statements\n')
  fprintf('  http://www.mathworks.com/help/techdoc/ref/try.html\n')
  fprintf('\n')
  fprintf(' MATLAB HELP\n')
  fprintf(' In the command window type: help help\n')
  fprintf(' In the command window type: help doc\n')
  fprintf(' In the command window type: help lookfor\n')
  fprintf('\n');
  fprintf('NOTE: Pressing Ctrl-C when Matlab is busy will throw an error\n');
  fprintf('and Matlab will stop after executing the current compiled binary\n');
  fprintf('operation. A compiled binary operation cannot be interrupted.\n');
  fprintf('Ctrl-C can be tricky because it is also the short cut for "Edit->Copy"');
  fprintf('and its functionality changes depending on whether or not Matlab');
  fprintf('is busy executing (and no text is selected in the command window).\n');
  keyboard
  
  % BREAKPOINTS FROM COMMAND LINE
  % =======================================================================
  % Breakpoints are useful when you want to debug a specific section of
  % code, but do not want to have to step through everything before that
  % section. This is especially important when there are many lines of code
  % being executed before the point where you want to debug such as large
  % for loops.
  
  % Clear all debug break points
  dbclear all;
  
  % Get the current function call stack point
  stack_info = dbstack
  
  % Set a breakpoint (can also be done by clicking just to the right of the
  % line number in the editor)
  dbstop('MATLAB_Tutorial',sprintf('%d',stack_info.line+10))
  
  % Clear a breakpoint (can also be done by clicking the red breakpoint
  % dot)
  dbclear('MATLAB_Tutorial',sprintf('%d',stack_info.line+10))
  
  x = 0;
  
  % BREAKPOINTS USING EDITOR
  % =======================================================================
  % Practice setting a breakpoint to skip a for loop
  stack_info = dbstack;
  dbstop('MATLAB_Tutorial',sprintf('%d',stack_info.line+6));
  % Instead of stepping through this for loop, run dbcont after placing the
  % breakpoint at "x = 0;"
  for x = 1:100000
    y = x^2;
  end
  x = 0;
  
  % CONDITIONAL BREAKPOINT
  % =======================================================================
  % Set a conditional breakpoint (can also be done by adding a breakpoint
  % and then right clicking on it and choosing "Set/Modify Condition...")
  stack_info = dbstack;
  dbstop('MATLAB_Tutorial',sprintf('%d',stack_info.line+7),'if','x==50000')
  dbstop('MATLAB_Tutorial',sprintf('%d',stack_info.line+9))
  % Instead of stepping through this for loop, run dbcont
  for x = 1:100000
    % When it stops at x == 50000, you can examine variables and then
    % run dbcont again to finish the for loop
    y = x^2;
  end
  x = 0;
  
  % TRY/CATCH
  % =======================================================================
  % try/catch lets you catch an error and deal with it instead of the
  % program execution stopping that occur and handle the error
  
  % This while loop shows a simple example of using try/catch
  while true
    try
      your_input = input('Enter a number less than 5 (enter 0 to exit loop): ');
      if your_input == 0
        break;
      end
      if your_input >= 5
        error('Your input was not less than 5.');
      else
        fprintf('Your input was less than 5.\n');
      end
    catch ME
      ME.getReport
    end
  end
  
  % For example, you can keep a bad line of code from breaking the program
  clear asdf;
  try
    asdf(1); % We know this will fail since this variable was just cleared
  end

  % DBSTOP IF ERROR
  % =======================================================================
  % A useful debugging tool is to have Matlab stop anytime an error happens
  % which is not caught by a try/catch block. Matlab does not allow the
  % program to continue executing, but it does put you into debug mode at
  % the broken line of code so you can inspect what went wrong.  Use
  % dbstack, dbup, and dbdown to navigate the function call stack.
  dbstop if error;
  
  % Try running the following after the error occurs in interp1
  % K>> dbstack
  % > In interp1 (line 161)
  %   In MATLAB_Tutorial (line 1917)
  %
  % K>> dbup
  % In workspace belonging to MATLAB_Tutorial (line 1917)
  %
  % K>> size(x)
  % ans =
  %      1     1
  
  x = 0; y = [1 2];
  interp1(x,y,2); % This fails because size(x) is too small
  
  fprintf('Section 8 Complete.\n');
  
  return
end

%% section 9: Signal processing
if section_number == 9
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 9: Signal processing\n');
  fprintf('========================================================\n')
  
  fprintf('Please review MATLAB_Tutorial_DataProcessing.ppt\n')
  keyboard
  
  % SINE WAVE IN NOISE
  % =======================================================================
  % fs: sampling frequency (Hz)
  fs = 200;
  % dt: sample interval (sec)
  dt = 1/fs;
  % Nt: number of time samples
  Nt = round(10/dt);
  % time: time axis (sec)
  time = dt * (0:Nt-1).';
  
  f0 = 13;
  phi = 20;
  x = 1.5 * cos(2*pi*f0*time + phi);
  
  figure(1); clf;
  plot(time, x);
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  xlim([0 0.6]);
  title('Sine wave without noise');
  
  % Add in noise (randn is Gaussian distributed noise with unity variance)
  x = x + 3*randn(Nt,1);
  
  figure(2); clf;
  plot(time, x);
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  xlim([0 0.6]);
  title('Sine wave with noise');
   
  
  % DISCRETE FOURIER TRANSFORM (DFT)
  % =======================================================================
  % xf: DFT of x
  xf = fft(x);

  % df: DFT frequency interval
  df = 1/(Nt*dt);
  % freq: DFT frequency axis
  freq = df*ifftshift(-floor(Nt/2):floor((Nt-1)/2));
  
  % fftshift, iffshift: shift zero frequency component to and from center
  % of vector.
  figure(3); clf;
  plot(fftshift(freq), fftshift(abs(xf)));
  grid on;
  xlabel('Frequency (Hz)');
  ylabel('DFT magnitude');
  title('Sine wave with noise');

  % Note the peak is at -13 and +13 Hz and is easily seen despite the noise
  
     
  % INTERPOLATION
  % =======================================================================
  % fs: sampling frequency (Hz)
  fs = 20;
  % dt: sample interval (sec)
  dt = 1/fs;
  % Nt: number of time samples
  Nt = round(1/dt);
  % time: time axis (sec)
  time = dt * (0:Nt-1).';
  
  f0 = 3;
  x = cos(2*pi*f0*time);

  % Without interpolation, the cosine is hard to see
  figure(1); clf;
  plot(time, x,'o');
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  title('Sine wave without noise');
  
  % interp1: linear interpolation
  % Mt: oversampling factor
  Mt = 10;
  dt_Mt = dt/Mt;
  Nt_Mt = Nt*Mt;
  time_Mt = dt_Mt * (0:Nt_Mt-1).';
  
  x_Mt = interp1(time,x,time_Mt,'linear','extrap');
  
  figure(1);
  hold on;
  plot(time_Mt, x_Mt,'r-');
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  title('Sine wave without noise');
  
  % interp1: spline interpolation
  % Mt: oversampling factor
  Mt = 20;
  dt_Mt = dt/Mt;
  Nt_Mt = Nt*Mt;
  time_Mt = dt_Mt * (0:Nt_Mt-1).';
  
  x_Mt = interp1(time,x,time_Mt,'spline','extrap');
  
  figure(1);
  hold on;
  plot(time_Mt, x_Mt,'g-');
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  title('Sine wave without noise');
  
  % interpft: Fourier domain interpolation (works only with Nyquist sampled
  % signals)
  x_Mt = interpft(x,Nt_Mt);
  
  figure(1);
  hold on;
  plot(time_Mt, x_Mt,'k-');
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  title('Sine wave without noise');
  
  ylim([-1.1 1.75]);
  legend('Sampled','Linear','Spline','Fourier','location','best');
  
  
  % FILTERING
  % =======================================================================
  % fs: sampling frequency (Hz)
  fs = 200;
  % dt: sample interval (sec)
  dt = 1/fs;
  % Nt: number of time samples
  Nt = round(10/dt);
  % time: time axis (sec)
  time = dt * (0:Nt-1).';
  
  f0 = 13;
  x = cos(2*pi*f0*time);
  
  figure(1); clf;
  plot(time, x);
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  xlim([2 2.6]);

  % Add in noise
  x = x + 3*randn(Nt,1);
  
  figure(1);
  hold on;
  plot(time, x,'r');
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  
  % Filter signal
  [B,A] = butter(2,(f0+[-1 1])/(fs/2));
  x = filtfilt(B,A,x);
  
  figure(1);
  plot(time, x,'g');
  grid on;
  xlabel('Time (sec)');
  ylabel('Voltage (V)');
  
  legend('cos','cos+noise','cos+noise+filt');

  
  % FITTING AND OUTLIERS
  % =======================================================================
  %polyfit, polyval, medfilt1, sgolayfilt
  
  
  % 2D INTERPOLATION
  % =======================================================================
  %delaunay

  
  % SMOOTHING IMAGES
  % =======================================================================
  %medfilt2
  %wiener2
  
  
  % WORKING WITH S-PARAMETER DATA (TRANSFER FUNCTIONS AND IMPULSE
  % RESPONSES)
  % =======================================================================
  % SXPParse, emquest_txt_reader, sparameters
  
end

%% section 10: Signal processing 2
if section_number == 10
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 10: Signal processing 2\n');
  fprintf('========================================================\n')
  
  keyboard
  
  % =======================================================================
  % Topic 1: Frequency Domain/Frequency Spectrum
  % =======================================================================
  % ALL signals can be represented as some combination of sine waves just
  % like the one you have created.  This means you can work in the frequency
  % domain to make some data processing easier.  In the frequency domain the
  % values along the x-axis represent sine waves of different frequencies
  % and the y-axis represents the power of each frequency component.  This
  % is called the frequency spectrum of the signal.  To view data in the
  % frequency domain in MatLab you can use the Fast Fourier Transform (fft)
  % function.
  
  % FFT
  % =======================================================================
  % The fft function runs more efficiently when the signal it is operating on
  % has a length that is a power of 2 (2^#).  The nextpow2 function returns
  % the next power of 2 length higher than the input length.
  Fs = 1000;                    % Sampling frequency
  T = 1/Fs;                     % Sampling interval
  % Create a time vector
  L = 100000;                     % Length of signal (number of samples)
  start_val = 0;                % Start at time 0
  stop_val = (L-1)*T;           % End after 1 second has elapsed
  t = linspace(start_val,stop_val,L); % Time vector
  f1 = 10; % Hz
  x3 = sin(2*pi*f1*t);
  NFFT = 2^nextpow2(L); % Input: 100000, Output: 17 --> 2^17 = 131072
  X = fft(x3,NFFT)/L;    % Fourier Transform of x with NFFT number of points
  freq = Fs*linspace(0,1,NFFT); % Frequency axis for X
  freq2= Fs/2*linspace(0,1,NFFT/2); % Frequency axis for Y
  
  % Plot double-sided amplitude spectrum
  figure(1); subplot(3,1,2); plot(freq2(1:5000),abs(X(1:5000)))
  title('Frequency Spectrum of a 10 Hz Sine Wave')
  xlabel('Frequency (Hz)')
  ylabel('|Y(f)|')
  
  
  % Notice the spike at 10 Hz, which is the frequency you used to create the
  % sine wave!
  
  % ifft
  % To reverse the process and return to the time domain you can use the
  % Inverse Fast Fourier Transform (ifft).
  x4 = ifft(X,NFFT)*L; % Return X to time-domain
  figure(1); subplot(3,1,3); plot(t(1:1000)*1e3,x4(1:1000))
  ylim([-1 1])
  title('10 Hz Sine Wave back in the Time Domain')
  xlabel('time (microseconds)')
  
  
  % The frequency domain signal has now been transformed back into a time
  % domain signal and is identical to our original signal!
  
  % Multiple Signals
  % Now add some more sine waves to the first one and see what happens to
  % the frequency spectrum of the signal.
  f2 = 20; f3 = 35;
  x5 = sin(2*pi*f1*t)+0.7*sin(2*pi*f2*t)+1.5*sin(2*pi*f3*t);
  
  figure(2); clf; subplot(2,1,1); plot(t(1:1000)*1e3,x5(1:1000))
  title('Signal Made of 10, 20, & 35 Hz Sine Waves')
  
  
  X2 = fft(x5,NFFT)/L;    % Fourier Transform of x with NFFT number of points
  
  % Plot double-sided amplitude spectrum
  figure(2); subplot(2,1,2); plot(freq2(1:5000),abs(X2(1:5000)))
  title('Frequency Spectrum of 10, 20, & 35 Hz Sine Waves')
  xlabel('Frequency (Hz)')
  ylabel('|Y(f)|')
  
  % fftshift/ifftshift
  % %
  % f3= Fs/2*linspace(-1,1,NFFT);  % Frequency axis for spectrum with zero-frequency at the center
  % figure(3); plot(f3,abs(fftshift(X)))
  % xlabel('Frequency (Hz)')
  % ylabel('|Y(f)|')
  % xlim([-200 200])
  
  % =======================================================================
  % Topic 3: Windowing/Filtering
  % =======================================================================
  % There are many instances when it is useful to apply a window (also called
  % a filter or taper) to data in order to remove some undesirable component.
  % Windows are particularly useful for reducing sidelobes or removing
  % unwanted noise from a signal.
  % In the time domain windows are applied using convolution.  Convolution is
  % applied by multiplying similarly delayed elements of a signal and filter
  % and summing them together.  In other words, if you have signal A which is
  % composed of {a1,a2,a3} and filter H which is composed of {h1,h2,h3} the
  % convolution of a and h (a*h) would be
  % {a1h1, a1h1+a2h2, a1h1+a2h2+a3h3, a2h2+a3h3, a3h3}.
  % In the frequency domain a filter can be applied by simply multiplying the
  % frequency spectrum of the signal by the frequency spectrum of the filter.
  
  % FILTERING
  % =======================================================================
  Fs2 = 120e6;                        % Sampling frequency
  t2 = 0:1/Fs2:10e-6;                 % Time vector
  NFFT2 = 2^nextpow2(length(t2));     % Number of points
  ch = chirp(t2,140e6,10e-6,160e6);   % Chirped signal
  
  %Add noise to the chirped frequency spectrum
  chi = ch + 0.2*sin(2*pi*155e6*t2);
  CHI = fft(chi); % Transform to frequcney domain
  f4 = linspace(120,180,600);     % Create frequency axis
  
  % Plot original spectrum
  figure(5); plot(f4,20*log10(abs(CHI(1:600)./CHI(300))))
  title('Chirped Pulse Frequency Spectrum')
  xlabel('Frequency (MHz)')
  ylabel('Spectrum Power (dB)')
  
  
  % Bandpass filter data from 145 to 155
  b = fir1(64,[0.40 0.60]);
  data = filter(b,1,chi);
  DATA = fft(data);
  
  figure(5); hold on;
  plot(f4,20*log10(abs(DATA(1:600)./DATA(300))),'r')
  
  
  % Bandpass filter data from 142 to 152
  b2 = fir1(64,[0.35 0.55]);
  data2 = filter(b2,1,chi);
  DATA2 = fft(data2);
  
  figure(5);
  plot(f4,20*log10(abs(DATA2(1:600)./DATA2(300))),'g')
  hold off
  
  
  % WINDOWING
  % =======================================================================
  % Rectangular Window
  L = 1000;
  win_len = 201;
  half_win = (win_len-1)/2;
  spectrm = zeros(1,L);
  spectrm(L/2-half_win:L/2+half_win) = 1;
  
  figure(3); clf; subplot(2,1,1); plot(spectrm)
  axis([300 700 -0.5 2])
  title('Frequency Domain')
  xlabel('Frequency (Hz)')
  
  Y3 = spectrm.*spectrm;
  % figure(3); subplot(3,1,2); plot(Y3)
  % ylim([-0.5 2])
  
  NFFT = 2^nextpow2(L);
  X3 = ifft(spectrm,NFFT)/L;
  figure(3); subplot(2,1,2); plot(lp(fftshift(X3/max(X3))))
  title('Time Domain')
  xlabel('Time (s)')
  axis([482 542 -30 5])
  
  % Triangular Window
  tri_win = zeros(1,L);
  tri_win(L/2-half_win:L/2+half_win) = triang(win_len);
  figure(3); subplot(2,1,1); hold on; plot(tri_win,'r')
  
  Y4 = spectrm.*tri_win;
  % figure(3); subplot(3,1,2); hold on; plot(Y4,'r')
  % ylim([-0.5 2])
  
  X4 = ifft(Y4,NFFT)/L;
  figure(3); subplot(2,1,2); hold on; plot(lp(fftshift(X4/max(X4))),'r')
  axis([482 542 -50 5])
  
  
  % Hanning Window
  hann_win = zeros(1,L);
  hann_win(L/2-half_win:L/2+half_win) = hanning(win_len);
  figure(3); subplot(2,1,1); hold on; plot(hann_win,'g')
  
  Y5 = spectrm.*hann_win;
  % figure(3); subplot(3,1,2); hold on; plot(Y5,'g')
  % ylim([-0.5 2])
  
  X5 = ifft(Y5,NFFT)/L;
  figure(3); subplot(2,1,2); hold on; plot(lp(fftshift(X5/max(X5))),'g')
  axis([482 542 -70 5])
  
  % Blackman Window
  blackman_win = zeros(1,L);
  blackman_win(L/2-half_win:L/2+half_win) = blackman(win_len);
  figure(3); subplot(2,1,1); hold on; plot(blackman_win,'c')
  
  Y5 = spectrm.*blackman_win;
  % figure(3); subplot(3,1,2); hold on; plot(Y5,'c')
  % ylim([-0.5 2])
  
  X5 = ifft(Y5,NFFT)/L;
  figure(3); subplot(2,1,2); hold on; plot(lp(fftshift(X5/max(X5))),'c')
  axis([482 542 -90 5])
  
  other_sig = 10^(-20/20)*blackman(20).^2;
  figure(3); subplot(2,1,2); hold on; plot([512:531],lp(other_sig),'k')
  
  
  % filter/filtfilt/fir1
  % =======================================================================
  Fs2 = 120e6;                        % Sampling frequency
  t2 = 0:1/Fs2:10e-6;                 % Time vector
  NFFT2 = 2^nextpow2(length(t2));     % Number of points
  ch = chirp(t2,140e6,10e-6,160e6);   % Chirped signal
  CH = fft(ch,NFFT2);                 % Transform to frequency domain
  CH = CH(1:NFFT2/2);
  len = NFFT2/2;                      % Length of the data
  uniform_win = ones(1,len);          % Remove extraneous data
  f3 = linspace(120,180,NFFT2/2);     % Create frequency axis
  % Uniform window
  tri_win     = triang(len).';        % Triangular window
  hann_win    = hann(len).';          % Hanning window
  
  % Apply uniform window to frequency-domain data and normalize
  chf_uni = uniform_win.*CH;
  chf_uni = chf_uni./max(chf_uni);
  
  figure(4);
  subplot(3,1,1); plot(f3,20*log10(abs(CH./max(CH)))); grid; axis tight
  xlabel('Frequency (MHz)')
  ylabel('Power (dB)')
  subplot(3,1,2); plot(uniform_win,'b'); grid;
  axis([0 len 0 1.5])
  subplot(3,1,3); plot(f3,20*log10(abs(chf_uni))); grid; axis tight
  xlabel('Frequency (MHz)')
  ylabel('Power (dB)')
  
  
  % Apply triangular window to frequency-domain data and normalize
  chf_tri = tri_win.*CH;
  chf_tri = chf_tri./max(chf_tri);
  
  subplot(3,1,2); hold on; plot(tri_win,'r')
  subplot(3,1,3); hold on; plot(f3,20*log10(abs(chf_tri)),'r'); grid; axis tight
  
  
  % Apply hanning window to frequency-domain data and normalize
  chf_hann = hann_win.*CH;
  chf_hann = chf_hann./max(chf_hann);
  
  subplot(3,1,2); hold on; plot(hann_win,'g')
  hold off
  subplot(3,1,3); hold on; plot(f3,20*log10(abs(chf_hann)),'g'); grid; axis tight
  hold off
  
  legend('Uniform','Triangle','Hanning')
  
  
  % =======================================================================
  % Topic 4: Other Useful Functions
  % =======================================================================
  
  % Interpolation
  % =======================================================================
  % Interpolation is method of constructing new data points within a set of
  % known data points.  For example, if you are paid at $10 per hour you know
  % that after one hour you have $10, two hours you have $20,... eight hours
  % you have $80.  If you worked 6 1/2 hours one day you would need to
  % interpolate to find the amount you are owed for that day because you
  % worked less than the sampling rate for the payscale (i.e. one hour).
  % Because the function of your wages over time (amount owed = 10*hours
  % worked) is a linear function you can interpolate between the values that
  % you know ($60 for 6 hours and $70 for 7 hours) to find the amount you are
  % owed for 6 1/2 hours worked which is $65.
  % Other reasons for using interpolation are to increase the definition of a
  % signal to see details that are currently hidden by the coarseness of the
  % resolution.  A note of warning, a signal has a finite resolution.
  % Interpolation can be used to reveal details up to the natural resolution
  % of a signal but cannot be used to improve the resolution.
  
  % Exercise 1:
  sig = linspace(1,8,8);
  time = 1:length(sig);
  figure(6); clf; plot(time,sig,'.')
  xlabel('Time (hours)')
  
  
  time_interp = .5:.5:length(sig);
  sig_interp = interp1(time,sig,time_interp,'linear','extrap');
  figure(6); hold on; plot(time_interp,sig_interp,'or')
  
  
  % Exercise 2:
  f0 = 100;
  fs = 1000; % Hz
  time3 = 0:1/fs:.1;
  L = length(time3);
  sig3 = sin(2*pi*f0*time3);
  n_s = sig3 + randn(1,length(time3));
  figure(7); clf; plot(time3,n_s,'-.')
  
  
  M = 2;
  highres_time = interp1(linspace(1,L,L),time3,linspace(1,L,M*L));
  highres_sig = interpft(n_s,M*L);
  figure(7); hold on; plot(highres_time,highres_sig,'-or')
  
  
  M = 20;
  higherres_time = interp1(linspace(1,L,L),time3,linspace(1,L,M*L));
  higherres_sig = interpft(n_s,M*L);
  figure(7); hold on; plot(higherres_time,higherres_sig,'-og')
  
  
  % Polyfit/Polyval
  % =======================================================================
  % Sometimes you will have a random signal that does not strcitly adhere to
  % any closed-form function.  In this case you may have to approximate the
  % closest polynomial.  You can use the polyfit function to determine the
  % polynomial coefficients that most closely fit the random signal.  Then
  % the polyval function can be used to recreate the function over the
  % desired range.
  
  t = -2:.01:2;
  rcof = randn(1,5); % Random polynomial coefficients
  poly_sig = rcof(1)*t.^4+rcof(2)*t.^3+rcof(3)*t.^2+rcof(4)*t-rcof(5);
  poly_noise = randn(1,length(t));
  sig = poly_sig+poly_noise;
  figure(8); clf; plot(t,sig)
  
  p1 = polyfit(t,sig,1);
  p2 = polyfit(t,sig,2);
  p3 = polyfit(t,sig,3);
  p4 = polyfit(t,sig,4);
  
  figure(8); hold on; plot(t,polyval(p1,t),'r')
  
  figure(8); hold on; plot(t,polyval(p2,t),'g')
  
  figure(8); hold on; plot(t,polyval(p3,t),'c')
  
  figure(8); hold on; plot(t,polyval(p4,t),'m')
  
  fprintf('Section 10 Complete.\n');
  return
  
end


%% Section 11: Signal processing 3: s-parameters
if section_number == 11
    % author: Audrey Evans
    
    
  fprintf('\n=========================================================\n');
  fprintf(' Section 11: Signal processing 3: s-parameters\n');
  fprintf('========================================================\n')
  
  keyboard
  
    % =======================================================================
    % S-Parameters
    % =======================================================================
    % Scattering parameters (s-parameters) are used to described electrical 
    % behavior of linear networks that are undergoing steady state stimuli.
    % These parameters are useful in electrical engineering, especially when
    % working with antenna, radio, and microwave engineering.

    close all; % close all existing figures

    % You can load s parameters two different ways; either using the matlab
    % toolbox or using SXPParse:

    % *** NOTE! Make sure that Accum_bpf.s2p is in your tmp folder.***

    % EXAMPLE 1: not using the matlab toolbox
    % load s parameters using SXPParse

    [freq,data]=SXPParse('C:\tmp\Accum_bpf.s2p');

    % EXAMPLE 2: using matlab toolbox
    % load s parameters using matlab toolbox

    s_obj = sparameters('C:\tmp\Accum_bpf.s2p');


    % MAGNITUDE VS. FREQUENCY OF S-PARAMETERS
    % =======================================================================

    figure(1);
    % First, you can plot the S11,S12,S21,and S22 parameters:
    subplot(2,2,1);
    plot(freq*1e-9,squeeze(db(data(1,1,:),'voltage'))); 
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    title('S_1_1 parameters');
    hold on;

    subplot(2,2,2);
    plot(freq*1e-9,squeeze(db(data(1,2,:),'voltage')));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    title('S_1_2 parameters');
    hold on;

    subplot(2,2,3); 
    plot(freq*1e-9,squeeze(db(data(2,1,:),'voltage')));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    title('S_2_1 parameters');
    hold on;

    subplot(2,2,4);
    plot(freq*1e-9,squeeze(db(data(2,2,:),'voltage')));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    title('S_2_2 parameters');
    hold off;

    figure(2)
    % It can be useful to plot S11 and S22 together, as well as S12 with S21:
    plot(freq*1e-9,squeeze(db(data(1,1,:),'voltage')));   hold on;
    plot(freq*1e-9,squeeze(db(data(2,2,:),'voltage'))); 
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    legend('S_1_1','S_2_2');
    title('S_1_1/S_2_2 parameters');

    figure(3)
    % S12 and S21 together:
    plot(freq*1e-9,squeeze(db(data(1,2,:),'voltage')));   hold on;
    plot(freq*1e-9,squeeze(db(data(2,1,:),'voltage'))); 
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    legend('S_1_2','S_2_1');
    title('S_1_2/S_2_1 parameters');


    % PHASE ANGLE VS. FREQUENCY OF S-PARAMETERS
    % =======================================================================

    S11 = squeeze(data(1,1,:));
    S12 = squeeze(data(2,1,:));
    S21 = squeeze(data(1,2,:));
    S22 = squeeze(data(2,2,:));


    figure(4)
    phaseangle = angle(S11); % This returns the phase angles in radians of the matrix.
    plot(freq*1e-9,phaseangle);  
    xlabel('Frequency (GHz)');
    ylabel('Phase (rad)');
    title('Phase angle plot for S_1_1');

    figure(5)
    phaseangle = angle(S21); % This returns the phase angles in radians of the matrix.
    plot(freq*1e-9,phaseangle);  
    xlabel('Frequency (GHz)');
    ylabel('Phase (rad)');
    title('Phase angle plot for S_2_1');


    % IMPULSE RESPONSE
    % =======================================================================
    % Next, we want to plot the impulse response of the s parameters.  An
    % impulse response is the plot of the behavior of the system when it is fed
    % a very brief input signal.  For this plot, we will need to convert the
    % frequency domain to the time domain.  In order to do this, we will use
    % ifft.

    figure(6)
    plot(db(S21,'voltage'))

    % Windowing for S21 in frequency domain
    H_idx = find(freq > 600e6 & freq < 900e6);
    H = hanning(length(H_idx));
    S21_window = zeros(size(S21));
    S21_window(H_idx) = S21(H_idx) .* H;

    df = freq(2)-freq(1);
    Nt = length(S21);
    BW = df*Nt;
    dt = 1/BW;
    time = dt*(0:Nt-1);

    % Convert from frequency domain to time domain:
    S21_td = ifft(S21_window);
    clf; % clear current figure
    plot(time*1e9, db(S21_td,'voltage'),'r'); hold on

    % TIME GATED S21 FREQUENCY PLOT
    % =======================================================================
    % source: http://incompliancemag.com/article/s-parameter-data-correction-using-time-domain-gating-for-pcb-and-cable-applications/
    % Time gating is a method used on impulse response functions that removes
    % reflections that are caused by end connectors or other discontinuities. 
    % We will plot the time gated S21 plot along with the impulse response that
    % we have already created.

    % Apply time domain window to S21_td
    H_idx = find(time > 0 & time < 20e-9);
    % We are using a hanning window.  See section 10 topic 3 of the matlab
    % tutorial.
    H = hanning(length(H_idx));  
    S21_td_window = zeros(size(S21_td));
    S21_td_window(H_idx) = S21_td(H_idx) .* H;
    plot(time*1e9, db(S21_td_window,'voltage'),'g');

    hold on
    S21_td = ifft(S21);
    plot(time*1e9, db(S21_td,'voltage'));
    xlabel('Time (ns)');
    ylabel('Magnitude (dB)');
    legend('S_2_1 = frequency domain window','S_2_1 = freq and time windowed','S_2_1 = time domain, not windowed');
    title('Impulse Response and Time Gating');
    % We will use axis() to zoom in so that we can see the impulse response
    % better.
    axis([-10,500,-160,-10]); 
    axis([-10,200,-160,-10]); 
    axis([-10,100,-130,-10]); 
    axis([-10,50,-110,-10]); 

    figure (7)
    % Now we want to compare the time gated plot, the windowed plot, and the
    % original plot
    S21_timegated = fft(S21_td_window);
    plot(time*1e9,db(S21_timegated,'voltage')); hold on
    plot(time*1e9,db(S21_window,'voltage')); hold on
    plot(time*1e9,db(S21,'voltage')); 
    xlabel('Time (ns)');
    ylabel('Magnitude (dB)');
    legend('S_2_1 time gated','S_2_1 frequency windowed','S_2_1');
    title('Impulse Response and Time Gating');



    % GROUP DELAY FOR S21
    % =======================================================================
    % source: http://www.microwaves101.com/encyclopedias/group-delay-measurements
    %   Group delay is the change in phase angle over the change in frequency.
    % In other words, it is the derivative of the phase angle with respect to
    % the frequency divided by 2pi.
    %   We are interested in the group delay for S21 because S21 is the
    % measurement of the complex output/input transfer function.  This means
    % that S21 represents both amplitude and phase between input and output
    % signals.
    %   First, we need to unwrap the phase angle using unwrap().  This changes
    % absolute jumps >= pi to their 2*pi compliment
    % We will then take the diff(unwrap(phaseangle))/diff(freq) vs freq.

    figure(8)

    x=unwrap(angle(data));

    % We have to take the transpose of the matrix here in order to be able to divide diff(phaseangle) by diff(freq)
    y=diff(freq).'; 
    plot(freq(1,1:20000)*1e-9, squeeze(diff(x(1,2,:))) ./ y .*1e4);

    % We will use axis() to zoom in on the figure because the function outside of this domain is just noise due to the unwrapping function not working.
    axis([0.3,1.1,-0.8,0.8]); 
    xlabel('Frequency (GHz)');
    ylabel('\partial\phi/\partialf x 1000/2\pi'); %Side note: using a '\pi' makes an actual pi symbol
    title('Group delay of S_2_1');




    % S-PARAMETER PLOTS USING THE MATLAB TOOLBOX
    % =======================================================================
    % source: http://www.mathworks.com/help/rf/ug/rfplot.html
    % rfplot will plot the magnitude in dB vs. frequency of all S-parameters,
    % while rfplot(s_obj,i,j) will plot Sij.

    figure(9)

    subplot(2,2,1)
    rfplot(s_obj,1,1); %plot S11
    xlabel('frequency');
    ylabel('mag. in dB'); hold on;

    subplot(2,2,2)
    rfplot(s_obj,1,2); %plot S12
    xlabel('frequency');
    ylabel('mag. in dB'); hold on;

    subplot(2,2,3)
    rfplot(s_obj,2,1); %plot S21
    xlabel('frequency');
    ylabel('mag. in dB'); hold on;

    subplot(2,2,4)
    rfplot(s_obj,2,2); %plot S22
    xlabel('frequency');
    ylabel('mag. in dB'); hold off;

    % You can also plot all four together by using rfplot(s_obj).

    figure(10)
    % S11 and S22 plotted together
    rfplot(s_obj,1,1,'r'); hold on;
    rfplot(s_obj,2,2,'b');


      fprintf('Section 11 Complete.\n');
  return
end


%% Section 12: Layer Plots
if section_number == 12
  fprintf('\n=========================================================\n');
  fprintf(' Section 12: Layer Plots \n');
  fprintf('========================================================\n')
  
  keyboard
  
  % SCATTER PLOT FOR OPS POLYGON
  % =======================================================================
  % This example shows how to query two layers from the Open Polar Server (ops)
  % database for a given region
  % See opsScatterPlotExample.m and opsGetPointswithinPolygon.m
  
  %These imputs are necessary for the opsGetPointsWithinPolygon function.
  %'rds' refers to the type of system used and 'artic' refers to the
  %location of the polygon region
  sys = 'rds';
  param = [];
  param.properties.location = 'arctic';
  
  % Polygons can be made at ops3.cresis.ku.edu (Select Artic and then
  % Draw Polygon. Close the polygon by double clicking then copy and
  % paste the WKT Polygon information into the form
  % param.properties='Polygon((...))'
  param.properties.bound='POLYGON((-32.92164206498864 68.61008511349851,-32.93658231364552 68.6167888873584,-32.95116092292497 68.61181999966429,-32.9254919803559 68.6078299519973,-32.92164206498864 68.61008511349851))';
  
  % The function opsGetPointsWithinPolygon will use WKT Polygon information to
  % obtain the points within the polygon. Note: the larger the polygon the
  % more time the function will use to obtain the information.
  fprintf('Getting points. This may take a minute (%s)\n', datestr(now));
  [status,data] = opsGetPointsWithinPolygon(sys,param);
  fprintf('  Done (%s)\n', datestr(now));
  
  % We need information about the area. By loading geotiff information, we
  % can use the 'geotiffinfo' command and 'projfwd' function to get more
  % information for the axes.
  geotiff_fn = ct_filename_gis([],fullfile('greenland','Landsat-7','Greenland_natural_150m.tif'));
  proj = geotiffinfo(geotiff_fn);
  [data.properties.X,data.properties.Y] = projfwd(proj,data.properties.Lat,data.properties.Lon);
  
  % Now we can make a scatterplot of the bed elevations
  figure (1); clf;
  hold on;
  scatter(data.properties.X/1e3, data.properties.Y/1e3, [], data.properties.Elevation - data.properties.Bottom);
  axes_elev = gca;
  h_colorbar=colorbar;
  set(get(h_colorbar,'YLabel'),'String','WGS-84 elevation(m)')
  axis normal;
  xlabel('X (km)')
  ylabel('Y (km)')
  
  % Alternately, we can plot the quality levels of the data
  proj = plot_geotiff(geotiff_fn,[],[], 2);
  hold on;
  
  % The for loop uses the quality levels (1, 2, or 3) for the bottum layer and
  % plots a data point associated with the level (green, yellow, or red).
  quality_color = {'g.' 'y.' 'r.'};
  for quality_level = [1 2 3]
    X = data.properties.X/1e3;
    Y = data.properties.Y/1e3;
    Q = data.properties.Bottom_Quality;
    quality_mask = Q == quality_level;
    if any(quality_mask)
      h(quality_level) = plot(X(quality_mask), Y(quality_mask), quality_color{quality_level});
    end
  end
  hold off;
  
  axes_quality = gca;
  axis normal;
  xlabel('X (km)');
  ylabel('Y (km)');
  legend('Good','Medium','Bad');
  linkaxes([axes_elev axes_quality],'xy');

  
  % LOADING LAYERS
  % =======================================================================
  % This example shows how to query layer data from ops for specific data
  % frames.
  % See opsLoadLayers.m and runOpsLoadlayers.m
  
  % We clear the structure param and specify relevant information (i.e.
  % season, radar name, date, etc.)
  param = [];
  param.season_name = '2011_Greenland_P3';
  param.radar_name = 'mcords2';
  param.day_seg = '20110414_04';
  param.cmd.frms = 1;
  param.post.ops.location = 'arctic';
  global gRadar;
  param = merge_structs(param,gRadar);
  
  % Before we load the data, we must the layers themselves (in this case
  % top and bottum) and their source (ops).
  layer_params = [];
  layer_params(1).name = 'surface';
  layer_params(1).source = 'ops';
  layer_params(2).name = 'bottom';
  layer_params(2).source = 'ops';
  
  % Now we can load the data from ops
  layers = opsLoadLayers(param,layer_params);
  
  % First, we plot the surface layer
  figure(1); clf;
  along_track = geodetic_to_along_track(layers(1).lat,layers(1).lon);
  h_plot = plot(along_track/10^3, layers(1).twtt*10^6);
  ylabel('Two Way Travel Time (microseconds)')
  xlabel('Along Track Distance (km)')
  title('Layers vs Time')
  
  % Now we plot the bottum layer
  hold on;
  along_track = geodetic_to_along_track(layers(2).lat,layers(2).lon);
  h_plot = plot(along_track/10^3, layers(2).twtt*10^6);
  hold off
  
  
  % ECHOGRAM WITH DEPTH AXIS
  % =======================================================================
  % This example shows how to load and plot an echogram with a depth axis
  % (rather than two way travel time) using a dielectric profile
  % See plot_L1B.m and elevation_compensation.m
  
  fn=fullfile(data_folder_path,'CRESIS_MAT_ECHODATA.mat');
  mdata=load_L1B(fn);
  
  % To correctly use the elevation_compensation function, we must clear
  % param and then specify certain inputs.
  param=[];
  
  % By specifying param.update is false (the default), we avoid making
  % changes to the mdata.Surface
  param.update_surf=false;
  
  % By specifying param.filter is true, we apply an along-track median
  % filter to the surface
  param.filter_surf=true;
  
  % This is the vector of dielectrics
  param.er_ice = linspace(1.5,3.5,101);
  
  % Next, we create a depth vector that will align with er_ice
  param.er_depth=linspace(0,100,101);
  
  % By specifying mode 2, we set the function to work towards a depth
  % axis rather than indices specified by param.depth
  param.elev_comp = 2;
  
  %Next, we specify the depth range
  param.depth = '[-200 1000]';
  
  % We call the function elevation_compensation which will give a depth
  % axis output
  [mdata_WGS84,depth_axis] = elevation_compensation(mdata,param);
  
  figure(2);
  clf;
  imagesc([],depth_axis,10*log10(mdata_WGS84.Data));
  xlabel('Range line');
  ylabel('Depth (m)')
  colormap(1-gray(256))
  
  
  %% Load echogram and layerData and sync + extract bottom power
  
  % Step 1: Load echogram file
  % echogram_fn: filename to L1B echogram data file, typical path might be:
  %   echogram_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Antarctica_DC8/CSARP_standard/20111014_07/Data_20111014_07_005.mat';
  echogram_fn = fullfile(data_folder_path,'CRESIS_MAT_ECHODATA.mat');
  % load_L1B: General cresis toolbox function that loads compressed
  % echograms and netcdf files in addition to the standard format.
  %mdata = load_L1B(fn);
  mdata = load(echogram_fn);
  
  % Step 2: Load layerData file
  % layer_fn: filename to layerData L2 layer file, typical path might be:
  %   layer_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Antarctica_DC8/CSARP_layerData/20111014_07/Data_20111014_07_005.mat';
  layer_fn = fullfile(data_folder_path,'CRESIS_MAT_DATA.mat');
  lay = load(layer_fn);
  
  % Step 3: Synchronize layer data with echogram data
  %  Units are two way travel time (twtt) in seconds
  %  Layer 1 is surface
  %  Layer 2 is bottom
  %  Value 2 is the full layer info, value 1 has just the manually entered points
  mdata.Surface = interp1(lay.GPS_time,lay.layerData{1}.value{2}.data,mdata.GPS_time);
  mdata.Bottom = interp1(lay.GPS_time,lay.layerData{2}.value{2}.data,mdata.GPS_time);
  
  % Step 4: Surface and Bottom converted to rows/range-bins from twtt
  mdata.Surface_Bin = interp1(mdata.Time,1:length(mdata.Time),mdata.Surface);
  mdata.Bottom_Bin = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);

  % Step 5: Plot echogram with synchronized layer data
  figure;
  imagesc(lp(mdata.Data));
  xlabel('Range lines');
  ylabel('Range bins');
  hcolor = colorbar;
  set(get(hcolor,'YLabel'),'string','Relative power (dB)');
  hold on;
  plot(mdata.Surface_Bin);
  plot(mdata.Bottom_Bin);
  hold off;
  
  % Step 6: Extract max bottom scattering power
  Nt = size(mdata.Data,1);
  Nx = size(mdata.Data,2);
  mdata.Bottom_Power_dB = zeros(1,Nx);
  % rbin_relative_search_rng: Search in a range around the bottom pick to
  % find the maximum power returned To return just the pixel power at the
  % bottom bin, set rbin_relative_search_rng = [0 0].
  rbin_relative_search_rng = [-5 10];
  for rline = 1:Nx
    rbin_search_rng = max(1,round(mdata.Bottom_Bin(rline))+rbin_relative_search_rng(1)) ...
      : min(Nt, round(mdata.Bottom_Bin(rline))+rbin_relative_search_rng(end));
    mdata.Bottom_Power_dB(rline) = max(10*log10(mdata.Data(rbin_search_rng, rline)));
  end

  % Step 7: Plot max bottom scattering power
  figure;
  plot(mdata.Bottom_Power_dB);
  xlabel('Range lines');
  ylabel('Bottom relative power (dB)');
  grid on;
  
  
end
    
    
