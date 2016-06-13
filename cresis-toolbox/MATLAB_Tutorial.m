% =========================================================================
%                      CReSIS MATLAB TUTORIAL
% =========================================================================

% 1. Follow steps here:
%    https://wiki.cresis.ku.edu/cresis/MATLAB_Tutorial
% 2. Run this tutorial by pressing F5 (choose "Change Folder" if asked)
%    Follow along in the terminal (command window) and editor

clc;
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
fprintf('9: data processing\n');

done = false;
while ~done
  try
    section_number = input('Enter tutorial to run [1-9]: ');
    section_number = round(section_number(1));
    if section_number >= 1 && section_number <= 9
      done = true;
    end
  end
end

if any(section_number == [4 5 6])
  global gdata_folder_path;
  if isempty(gdata_folder_path)
    gdata_folder_path = 'C:\tmp\Tutorial\Data\';
  end
  fprintf('Default Data Path: %s\n', gdata_folder_path);
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

clc
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
  % ---------------------------------------------------------
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
  % ---------------------------------------------------------
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
  % ---------------------------------------------------------
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
  % ---------------------------------------------------------
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
  % ---------------------------------------------------------
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
  
  fprintf(' Below we will give 3 file reading examples.\n')
  fprintf(' 1) Read a standard CReSIS CSV\n')
  fprintf(' 2) Read a CReSIS MAT file\n')
  fprintf(' 3) Read a binary lidar file (LVIS Lidar)\n')
  fprintf('\n')
  fprintf(' When you donwloaded the tutorials ZIP there was a "Data" folder\n')
  fprintf(' included. You will now be using that data.\n')
  fprintf('\n')
  fprintf(' IMPORTANT !!!!\n')
  fprintf('\n')
  
  keyboard
  
  
  fprintf('READ A STANDARD CReSIS CSV Using TEXTSCAN:');
  
  % Set the file path
  file_path = fullfile(data_folder_path,'CRESIS_CSV_DATA.csv');
  
  % Open the file in the editor to view the contents
  %  - This file is a "comma separated variable" or csv file
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
  %   o headerlines: there is 1 line that it is header
  csv_data = textscan(fid,'%f%f%f%f%f%f%f%f%s','delimiter',',','headerlines',1);
  
  % Close the data pointer (DONT FORGET THIS)
  fclose(fid);
  
  % Note that csv_data is a cell array
  % csv_data{1} is the first column of data from the file
  csv_data
  
  % This file contains ice thickness data from a flight track
  figure; % Create a new figure
  plot(csv_data{2}, csv_data{1}); % plot the longitude (2nd column) and latitude (1st column)
  xlabel('Longitude (deg)'); % add an xlabel to the plot
  ylabel('Latitude (deg)'); % add a ylabel to the plot
  
  keyboard;
  
  fprintf('READ A CReSIS MAT file:');
  
  % Set the file path
  file_path = fullfile(data_folder_path,'CRESIS_MAT_DATA.mat');
  mat_data = load(file_path);
  
  % Note this loads a standard CReSIS file called Layer Data. The format is
  % explained here:
  % https://wiki.cresis.ku.edu/cresis/CReSIS_Data_Format/s
  
  keyboard;
  
%   fprintf('READ A binary Lidar File (LVIS Lidar):');
%   
%   % WARNING: THE LOAD ON THIS PART MAY TAKE A SECOND. IT IS A LARGE BINARY
%   % FILE.
%   
%   % Set the file path
%   file_path = fullfile(data_folder_path,'LVIS_LIDAR_DATA.lge');
%   
%   % Specify the binary format information (This would be given in a README)
%   recordType = {'ulong','ulong','double','double','double','float','float','float','float','float'};
%   recordLen = [4 4 8 8 8 4 4 4 4 4];
%   
%   % Create a holder for the data to read in (see help numel())
%   bin_data = cell(1,numel(recordType));
%   
%   % Open a binary file pointer. Noter the 'rb' (See help fopen)
%   fid = fopen(file_path,'rb');
%   
%   % Loop through and use fseek and fread (See help on both) to read the
%   % data.
%   for b_idx = 1:numel(recordType);
%     fseek(fid, sum(recordLen(1:b_idx-1)),'bof');
%     bin_data{b_idx} = fread(fid, Inf, [recordType{b_idx}], sum(recordLen)-recordLen(b_idx),'ieee-be');
%   end
%   
%   % Close the data pointer (DONT FORGET THIS)
%   fclose(fid);
%   
%   % Now explore bin_data which is a cell of the binary data.
  
  keyboard;
  fprintf('Lets go through a simple example....');
  
  % Lets load some example data from a TXT file.
  
  % The file looks like this:
  
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
  
  % What if we want to know all the different grades the class has recieved?
  % Did everyone get A's or some B's and C's?
  
  fprintf('Lets use the unique function ....');
  unique_grades = unique(grade_data{3});   % NOTE THE 3rd cell {3} is the list of grades from the file. Make sense?
  
  % If we look at "unique_grades" we can see that the grades of "A,B,C,and F"
  % were recieved, but not D's
  
  % Next lets get the names of the students that recieved "F" using logical
  % indexing.
  
  f_indexes = logical(strcmpi(grade_data{3},'F'));    % STRCMPI compares a string and ignores case
  
  % A "logical" array is 1's and 0's. A "1" represents a true condition a "0"
  % a false condition. The above will give us an array of 1's and 0's the
  % length of the grades array with 1's at the indexes of the f's. We can use
  % this index to get the name/s of the students with "f" for a grade.
  
  f_students = grade_data{1}{f_indexes};      % Keep only names from grade_data{1} that have "1" in the "f_indexes"
  
  % Looks like only "Joe" got an F!
  
  % Say some professor wants a list of all of the "A" students in a new TXT
  % file. We now know how to do this!
  
  a_indexes = logical(strcmpi(grade_data{3},'a'));    % Find all the "A" or "a" indexes.
  
  % Now lets use the indexes to get the names and I'ds of all "A" students.
  % We can use a loop to do this for every column in the cell "grade_data"
  
  for data_idx = 1:length(grade_data)
    grade_data{data_idx} = grade_data{data_idx}(a_indexes);
  end
  
  % "grade_data" now only contains the information for the A students.
  
  keyboard;
  fprintf('  Lets write this to a new file ....');
  
  new_grade_file = fullfile(data_folder_path,'A_GRADES_TXT_DATA.txt');
  
  fid = fopen(new_grade_file,'w+');   % The 'w+' means fopen will create a new file or discard existing contents of a file.
  
  % First lets print the header.
  fprintf(fid,'%s,%s,%s\n','NAME','USERID','GRADE');
  
  % Now we can use a loop to print the data one line at a time.
  for data_idx = 1:length(grade_data{1})
    fprintf(fid,'%s,%d,%s\n',grade_data{1}{data_idx},grade_data{2}(data_idx),grade_data{3}{data_idx});
  end
  
  fclose(fid);
  
  % SAVE
  % ------------------------------------
  
  % Save is a very simple command that saves various variables and
  % workspaces in mat or ascii format. See help save for more information.
  
  % That's it for section 4, you should now understand how to read and
  % write various file formats. Dont forget you can always use the help for
  % other file types!
  
  
  % That's it! Write some statements, make some loops and get comfortable!
  % Conditionals and loops are the gateway to programming!
  
  fprintf('Section 4 Complete.\n');
  return
end

%% Section 5: visualization
if section_number == 5
  
  fprintf('\n=========================================================\n');
  fprintf(' Section 5: visaulization\n');
  fprintf('========================================================\n')
  
  fprintf(' Here are a few videos to watch before this section:\n')
  fprintf(' http://www.mathworks.com/videos/using-basic-plotting-functions-69018.?s_tihtmld=srchtitle\n')
  fprintf(' http://www.mathworks.com/videos/visualizing-data-with-matlab-68917.html\n')
  
  fprintf(' Figures and Handles\n')
  fprintf(' ---------------------------------------------------------\n')
  
  fprintf(' Note the figure and handle are the basic manipulators for plots and\n')
  fprintf(' figures. Read through the documentation below for detailed\n')
  fprintf(' information...\n')
  
  fprintf(' http://www.mathworks.com/help/techdoc/learn_matlab/f3-15974.html\n')
  
  fprintf(' Lets create a few figures...\n')
  close all
  keyboard
  
  figure(1)
  % Creates a new figure with the handle 1
  
  h = figure('Name','MyPlot','NumberTitle','off');
  % Creates a figure with a custom title and stores the figure in handle h.
  
  clear h; close all;
  
  % Plot and Stem
  fprintf('\n') % ---------------------------------------------------------
  
  x = 1:.01:1.5;
  y = cos(2*pi*x);
  
  figure(1);      % creates a new figure
  plot(x,y);      % plot a cosine curve on figure 1
  
  
  figure(2);      % make another figure
  stem(x,y);      % stem a cosine curve on figure 1
  
  % Plot a CReSIS elevation on a time axis
  file_path = fullfile(data_folder_path,'CRESIS_MAT_DATA.mat');
  mat_data = load(file_path);
  
  figure(3);
  plot(mat_data.GPS_time,mat_data.Elevation,'r.');
  grid on;
  
  % Note the 'r.' This means plot point with dots colored red.
  % What if we change it to 'g-' ???
  
  figure(4)
  plot(mat_data.GPS_time,mat_data.Elevation,'g-');
  
  % See help plot for all the variations of color and symbol.
  
  clear x y mat_data; close all;
  
  
  % Images
  % ---------------------------------------------------------
  
  % We can plot CReSIS data as an image to view the "echogram"
  file_path = fullfile(data_folder_path,'CRESIS_MAT_ECHODATA.mat');
  echo_data = load(file_path);
  
  % Plot the image using a base 10 log plot db()
  figure(1)
  imagesc(db(echo_data.Data,'power'))
  
  % Plot using default x-axis (column number) and two way travel time for
  % the y-axis
  figure(2)
  imagesc([],echo_data.Time*1e6,db(echo_data.Data,'power'))
  hcolor = colorbar;
  set(get(hcolor,'YLabel'),'String','Relative power (dB)');
  xlabel('Range lines');
  ylabel('Two way travel time (us)');
  
  % Change the colormap
  colormap(gray(256)); % 256 color gray scale
  colormap(1-gray(256)); % invert gray scale
  colormap(hsv(256)); % hue only colormap (used to plot phase/angle)
  colormap(jet(256)); % "jet" colormap (red is big, blue is small)
  
  clear file_path echo_data; close all;
  
  
  % Sub-Plots
  % ---------------------------------------------------------
  
  % Sub plots allow us to use a single figure window to plot multiple
  % objects on different axis.
  
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
  
  clear file_path echo_data; close all;
  
  % Advanced Plots (More than one feature)
  % ---------------------------------------------------------
  
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
  % This is a faster way to do the same thing.
  [lat,lon,utc_time,thickness,elev,frame,surface,bottom,quality] = deal(csv_data);
  
  % Now we can plot surface and bottom. It's important to note in CReSIS data
  % variables surface and bottom are the distance from elevation to each
  % respective surface. We need to calculate the actual values.  
  a_surf = csv_data{5}-csv_data{7};
  a_bed = csv_data{5}-csv_data{8};

  figure(1);
  
  % We now want to "hold" the figure
  hold on;
  
  
  plot(csv_data{3},a_surf,'r-')
  plot(csv_data{3},a_bed,'g-')
  
  % Lets put a legend,title,and labels on the plot.
  title('Surface and Bottom Layers');
  xlabel('UTC Time (SOD)');
  ylabel('Elevation WGS84 (m)');
  legend('Surface','Bottom');   % Go in the order of the plots.
  
  clear file_path fid csv_data a_surf a_bed; close all;
  
  % Statistical Plots
  % ---------------------------------------------------------
  
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
  
  
  figure(3);
  x = 1:6;
  y = ceil(rand(1,10000)*6);
  hist(y,x);
  
  clear x y; close all;
  
  % Figure managment and manipulation
  % ---------------------------------------------------------
  
  % There are a set of commands that are used to acess and manipulate
  % figures through figure handles. They are:
  % get, set, gcf, gca, clf, cla, close.
  
  fprintf('Lets re-create the CReSIS sub-plot figure');
  file_path = fullfile(data_folder_path,'CRESIS_MAT_ECHODATA.mat');
  echo_data = load(file_path);
  figure(1)
  sp1 = subplot(2,2,1); imagesc(echo_data.GPS_time,echo_data.Depth,db(echo_data.Data,'power')); title('Echogram');
  sp2 = subplot(2,2,2); plot(echo_data.GPS_time); title('GPS Time')
  sp3 = subplot(2,2,3); plot(echo_data.GPS_time,echo_data.Surface,'k-'); title('Surface')
  sp4 = subplot(2,2,4); plot(echo_data.GPS_time,echo_data.Bottom,'k-'); title('Bottom')
  
  fprintf('Lets "get" the properties of the current axes');
  
  get(gca)  % This displays all of the properties of the current axes.
  
  % Not the current axes is that of the last plotted object, in this case
  % sp4 (The "Bottom" plot)
  
  fprintf('Lets turn the X and Y grid on using set');
  
  set(gca,'YGrid','on','XGrid','on')  % You could do this for all of the properties you saw in "get"
  
  fprintf('Lets set the y-limits of the echogram to [0 4000]');
  
  set(sp1,'YLim',[0 4000])    % Notice we pass the exact figure handle (sp1) intead of gca.
  
  fprintf('Lets turn the X and Y grid on for the surface plot');
  
  set(sp3,'YGrid','on','XGrid','on')
  
  fprintf('Lets clear only the GPS time axis');
  
  cla(sp2)
  
  fprintf('Finally, lets clear the entire figure, pause and then close the figure.');
  
  clf;
  pause(3); % Pause for 3 sec (So you can see the figure was cleared)
  close;
  
  
  fprintf('Color Manipulation\n')
  fprintf('---------------------------------------------------------\n')
  
  fprintf('We will now load the echogram plot and do some color manipulation...');
  
  figure(1)
  echo_f = imagesc(echo_data.GPS_time,echo_data.Depth,db(echo_data.Data,'power')); title('Echogram');
  echo_axis = axis; % Get the current axis properties
  axis([echo_axis(1) echo_axis(2) 0 4000]); % Set the new y-axis
  
  fprintf('Lets add a colorbar');
  
  colorbar;
  
  % Note the rest of the echogram is still part of the plot (the bottom
  % part) we just have limited our axes so you can see it.

  % To see the axis limits of the colorbar, we run caxis:
  caxis
  
  fprintf('Lets modify the color axes.');
  caxis([-120 -20]);
  caxis
  
  fprintf('Experiment with colormap now....')
  
  doc colormap;
  %   colormap (gray);
  %   colormap (1-gray)
  
  
  close;
  % Interactive Plots
  % ---------------------------------------------------------
  
  % The ginput() function allows a user to draw input to a figure which
  % matlab will store in a variable.
  
  % The following is a simple example that will allow a user to draw 5
  % points on a new figure.
  
  fprintf('Draw 5 points by clicking on the figure.');
  
  figure(1);                % Create a new figure
  [X Y] = ginput(5);        % Get input 5 times, store in X,Y
  plot(X,Y,'r+')            % Plot the input with a red +
  
  % Saving Plots
  % ---------------------------------------------------------
  
  % See the following for specifics on saving images in MATLAB
  % doc print
  % doc saveas
  
  
  close;
  fprintf('Lets load a built-in MATLAB image...');
  
  C = load('clown');    % Loads the MATLAB file
  Cl = image(C.X);
  
  fprintf('Lets save the clown to your DATA folder:');
  
  out_file = fullfile(data_folder_path,'clown_image')
  saveas(Cl,out_file,'jpeg');
  
  
  % 3D Plotting
  % ---------------------------------------------------------
  % See doc plot3 for more information.
  
  fprintf('Plot a 3D Helix');
  
  
  figure(1);              % Plot a helix. With 3D plots, the MATLAB figure
  t = 0:pi/50:10*pi;      % window tools are very useful: use the Hand tool
  plot3(sin(t),cos(t),t)  % to move the plot around, the Rotate tool to view it
  axis square;            % from different angles, and the Zoom tool to adjust
  grid on                 % the plot scale. They are located right below the menu
  % toolbar.
  
  
  % Surface Plots
  % ---------------------------------------------------------
  % See doc mesh; doc surf; doc contour;
  
  fprintf('Plot a 3D Mesh and 2D contours');
  
  figure(1);
  [X,Y] = meshgrid(-3:.125:3);    % Create evenly-spaced X and Y grids
  Z = peaks(X,Y);                 % Generate a Z grid using the peaks function
  meshc(X,Y,Z)                    % Create a mesh plot with a contour map below
  axis([-3 3 -3 3 -10 5])
  
  fprintf('Plot a 3D Mesh and 2D curtain');
  
  [X,Y] = meshgrid(-3:.125:3);
  Z = peaks(X,Y);
  meshz(X,Y,Z)                    % Same as above, but now a curtain is added instead of a contour map
  
  colormap hsv                    % Set hsv colormap
  axis([-3 3 -3 3 -10 5])
  
  fprintf('Plot a 2D peaks contour with labels');
  
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
  
  text_handle = clabel(C,h);      % Turn on contour labels so that the levels
  % of each contour are described
  % quantitatively
  
  set(text_handle,'BackgroundColor',[1 1 .6],'Edgecolor', [.7 .7 .7])
  
  
  % That's it! Plot away!
  
  fprintf('Section 5 Complete.\n');
  return
end

%% Section 6: mapping toolbox
if section_number == 6
  close all
  
   fprintf('\n=========================================================\n');
  fprintf(' Section 6: mapping toolbox\n');
  fprintf('========================================================\n')
  
  keyboard
  
  fprintf('Load and Plot World Coastline Data');
  
  load coast;         % Loads built in coastline vector data
  
  axesm mercator;     % Sets the projection to mercator
  framem;             % Adds a frame to the map window
  plotm(lat,long);    % Plots the lat/long of the coastline data.
  
  
  fprintf(' Clean Up:\n')
  clear lat long; close all;
  
  
  
  fprintf('Load the Mapping GUI (See code comments)');
  
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
  % 5) Shows the coasline vector data
  
  % If you cannot get your map to work correctly, or you are not sure what it
  % should look like use the code below to generate the map, then close and
  % re-try creation with the GUI.
  
  
  fprintf('Load the correct map w/o the GUI');
  
  load coast;
  axesm('MapProjection','mercator','Frame','on','Grid','on',...
    'MLineLocation',60,'PLineLocation',30,'MeridianLabel','on',...
    'ParallelLabel','on','MLabelLocation',60,'PLabelLocation',30);
  plotm(lat,long);
  
  % Clean Up:
  clear lat long; close all;
  
  
  fprintf('Load the korea grid and plot a texture map');
  
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
  
  
  fprintf('Load and plot a map with composite data types');
  
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
  
  
  fprintf('Load and plot a 3D map');
  
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
  
  
  fprintf('Read and explore and ESRI Shapefile');
  
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
  
  
  fprintf('Read and explore a GeoTIFF');
  
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
  
  
  fprintf('Load a GeoTIFF and Shapefile, plot the composite.');
  
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
  
  
  
  fprintf('Lets write some XY data to an ESRI shapefile...');
  
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
  
  keyboard;
  % Clean Up
  clear ans ax k_bbox ksmo_bound m_bbox mapstruct out_shp xk xm yk ym; close all;
  
  
  
  fprintf('Lets do some coordinate conversion ...');
  
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
  
  
  fprintf('Lets use the geodetic_to_along_track function ...');
  
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
  
  keyboard;
  % Clean Up
  clear along_track_dist coastline dist_between_coords elev_ft elev_m inspect_lat inspect_lon lat lon mstruct x y; close all;
  
  
  
  fprintf('Lets do a coordinate projection ...');
  
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
  
  % Lets project into the Transverse Mercator projection.
  proj = defaultm('tranmerc');
  [x,y] = projfwd(proj,lat,lon);
  
  % Great now we have the follwoing coordinate pair:
  % x = -1.0312 , y = 1.6838
  
  % What is we started with the Transverse Mercator coordinates and wanted
  % LAT/LON in degrees, WGS1984?
  
  [lat2,lon2] = projinv(proj,x,y);
  
  % So projinv is simply the inverse of projfwd. You can confirm that lat2 =
  % lat and lon2 = lon if you dont believe it!
  
  keyboard;
  % Clean Up
  clear lat lon x y lat2 lon2 proj;
  
  
  
  fprintf('Lets do a coordinate transformation ...');
  
  
  % Coordinate Transformation (mfwdtran,minvtran)
  % -------------------------------------------------------------------------
  % Say we have th following WGS84 Coordinate pair.
  
  lat = 38.9521;
  lon = -95.2640;
  
  % Lets transform these cordinates into UTM map coordinates. First we need
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
  
  keyboard;
  % Clean Up
  clear lat lon x y lat2 lon2 mstruct;
  
  
  fprintf('Lets apply all of the above to CReSIS data ...');
  
  % NOTE THIS EXAMPLE REQUIRES YOU TO BE ON A WINDOW MACHINE OR CHANGE TO
  % UNIX PATHS. ASK IF YOU DONT UNDERSTAND THIS.
  
  % DATA PATHS (MODIFY IF NECCESARY)
  
  csv_path = 'Z:\mdce\mcords\2009_Antarctica_DC8\CSARP_post\csv\Browse_2009_Antarctica_DC8.csv';
  gtiff_path = 'P:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_peninsula.tif';
  
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
  
  keyboard;
  
  % That's it! Map away!
  
  fprintf('Section 6 Complete.\n');
  return
end

%% Section 7: vectorization
if section_number == 7
  close all
  
   fprintf('\n=========================================================\n');
  fprintf(' Section 7: vectorization \n');
  fprintf('========================================================\n')
  
  fprintf(' See the following video and page for detailed information about improving\n')
  fprintf(' performance in MATLAB. This is important as you begin to work with\n')
  fprintf(' large datasets.\n')
  fprintf('Read this\n')
  fprintf('http://www.mathworks.com/company/newsletters/articles/programming-patterns-maximizing-code-performance-by-optimizing-memory-access.html\n')


  
  fprintf('This section will explain some methods of optimization using MATLAB.\n')
keyboard

fprintf('Now run dbcont to exit debug mode.');
  
  fprintf('Lets run a for loop with no preallocation:');
  
  clear x1 x2;
  len  = 10e6;
  tic
  for idx0 = 1:len
    x1(idx0) = cos(2*pi*100*idx0); % Data is created and stored on the fly.
  end
  toc
  
  
  fprintf('Lets run a for loop with preallocation:');
  
  tic
  x2 = zeros(1,len);    % We have pre-allocated a set of zeros to hold the data.
  for idx0 = 1:len
    x2(idx0) = cos(2*pi*100*idx0);
  end
  toc
  
  
  
  
  fprintf('Lets get rid of the loop and use a vectorized function:');
  
  tic
  t = 1:len;
  x3 = cos(2*pi*100*t);   % Here we use the cos() function. This function is vectorized.
  toc
  keyboard
  
  fprintf('Lets do the same thing, but with a MATRIX:');
  
  len = 100;
  
  r1 = randn(len);
  r2 = randn(len);
  fprintf('Now run dbcont again.');
  
  
  
  fprintf('For loop without preallocation')
  tic
  for idx0 = 1:len
    for idx1 = 1:len
      x1(idx0,idx1) = r1(idx0,idx1)*r2(idx0,idx1);
    end
  end
  toc
  
  fprintf('For loop with preallocation')
  tic
  x2 = zeros(len);
  for idx0 = 1:len
    for idx1 = 1:len
      x2(idx0,idx1) = r1(idx0,idx1)*r2(idx0,idx1);
    end
  end
  toc
  
  keyboard
  
  % Matrix multiply example
  tic
  x3 = r1.*r2;
  toc
  
  fprintf('Section 7 Complete.\n');
  return
end

%% Section 8: debugging
if section_number == 8
  
   fprintf('\n=========================================================\n');
  fprintf(' Section 8: debugging\n');
  fprintf('========================================================\n')
  
  fprintf('\n')
  fprintf('\n')
  fprintf(' There is no code here, rather a series of videos will introduce you to\n')
  fprintf(' debugging code in MATALB.\n')
  fprintf('\n')
  fprintf(' Here is another way of looking at "Debugging"\n')
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
  keyboard
  
  % BREAKPOINTS
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
  
  % Practice setting a breakpoint to skip a for loop
  stack_info = dbstack;
  dbstop('MATLAB_Tutorial',sprintf('%d',stack_info.line+6));
  % Instead of stepping through this for loop, run dbcont
  for x = 1:100000
    y = x^2;
  end
  x = 0;
  
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
  fprintf(' Section 9: signal processing\n');
  fprintf('========================================================\n')
  
  fprintf('look at the source code for more information\n')
  keyboard
  
  % NOTE THERE IS A SECOND POWERPOINT THAT GOES ALONG WITH THIS TUTORIAL,
  % REVIEW THE Data_Processing_Tutorial.pptx now.
  
  % This section is designed to introduce you to the some of the concepts
  % and frequently used functions associated with data processing at CReSIS.
  % We will not be delving into the theory or math behind these function,
  % just common practices and examples.
  
  %This is some boiler plate code to prepare the MatLab environment for
  %efficient use.
  format compact; format short
  
  % Topic 1: For Loops vs Vectorization (and preallocation of memory
  % There are mutliple ways to solve most coding problems but they are not
  % all created equal.  One of the most common issues you will face when
  % using MatLab is when to use For Loops and when to use vectorization.  The
  % main considerations with this question are speed, memory usage, and
  % computational complexity.
  % FOR LOOPS
  %   Advantages:
  %      -Requires less memory because only one element/row is being operated
  %       on at a time.
  %      -Only operates on one element/row at a time, which sometimes can
  %       is required based on the type of calulations you are performing.
  %   Disadvantages:
  %      -Significantly slower because it is only operating on one
  %       element/row at a time.
  %
  % VECTORIZATION
  %   Advantages:
  %      -Much faster because multiple elements/rows are being operated on
  %       simultaneously.
  %   Disadvantages:
  %      -Uses more memory because multiple elements/rows are being operated
  %       on simultaneously.
  %
  % Your choice will likely come down to your speed requirements or your
  % memory constrains.
  %
  % MEMORY PREALLOCATION
  % for and while loops that incrementally increase, or grow, the size of a
  % data structure each time through the loop can adversely affect
  % performance and memory use. Repeatedly resizing arrays often requires
  % that MATLAB spend extra time looking for larger contiguous blocks of
  % memory and then moving the array into those blocks. You can often improve
  % on code execution time by preallocating the maximum amount of space that
  % would be required for the array ahead of time.
  
  % Exercise 1
  % Goal: Create a sine wave using a for loop and vectorization and compare
  %       the amount of time need to create it using each method.
  % Functions used:
  %    linspace(start_val, stop_val, L)
  %       Creates a linearly spaced vector of length L which goes from
  %       start_val to stop_val.
  %    tic/toc: tic starts a counter and toc ends the counter and outputs the
  %             amount of time that has elapsed between tic & toc.
  Fs = 1000;                    % Sampling frequency
  T = 1/Fs;                     % Sampling interval
  % Create a time vector
  L = 100000;                     % Length of signal (number of samples)
  start_val = 0;                % Start at time 0
  stop_val = (L-1)*T;           % End after 1 second has elapsed
  t = linspace(start_val,stop_val,L); % Time vector
  
  % Create a 50 Hz sine wave
  f1 = 10; % Hz
  
  % Using a For loop w/o preallocation
  tic
  for idx = 1:L
    x1(idx) = sin(2*pi*f1*t(idx));
  end
  toc
  
  % Using a For loop with preallocation
  tic
  x2 = zeros(1,L);
  for idx = 1:L
    x2(idx) = sin(2*pi*f1*t(idx));
  end
  toc
  
  % Using vectorization
  tic
  x3 = sin(2*pi*f1*t);
  toc
  
  % Notice that the process is over twice as fast using preallocation with a
  % for loop and almost seven times faster using vectorization!
  
  figure(1); clf; subplot(3,1,1); plot(t(1:1000)*1e3,x3(1:1000))
  title('10 Hz Sine Wave')
  xlabel('time (microseconds)')
  
  
  
  % Topic 2: Frequency Domain/Frequency Spectrum
  % ALL signals can be represented as some combination of sine waves just
  % like the one you have created.  This means you can work in the frequency
  % domain to make some data processing easier.  In the frequency domain the
  % values along the x-axis represent sine waves of different frequencies
  % and the y-axis represents the power of each frequency component.  This
  % is called the frequency spectrum of the signal.  To view data in the
  % frequency domain in MatLab you can use the Fast Fourier Transform (fft)
  % function.
  
  % Exercise 2: fft
  % The fft function runs more efficiently when the signal it is operating on
  % has a length that is a power of 2 (2^#).  The nextpow2 function returns
  % the next power of 2 length higher than the input length.
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
  
  % Topic 3: Windowing/Filtering
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
  
  % Exercise 3: Filtering
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
  
  
  % Exercise 4:
  %
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
  
  
  % Blackman WindowWindow
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
  % Fs2 = 120e6;                        % Sampling frequency
  % t2 = 0:1/Fs2:10e-6;                 % Time vector
  % NFFT2 = 2^nextpow2(length(t2));     % Number of points
  % ch = chirp(t2,140e6,10e-6,160e6);   % Chirped signal
  % CH = fft(ch,NFFT2);                 % Transform to frequency domain
  % CH = CH(1:NFFT2/2);
  % len = NFFT2/2;                      % Length of the data
  % uniform_win = ones(1,len);          % Remove extraneous data
  % f3 = linspace(120,180,NFFT2/2);     % Create frequency axis
  %   % Uniform window
  % tri_win     = triang(len).';        % Triangular window
  % hann_win    = hann(len).';          % Hanning window
  %
  % % Apply uniform window to frequency-domain data and normalize
  % chf_uni = uniform_win.*CH;
  % chf_uni = chf_uni./max(chf_uni);
  %
  % figure(4);
  % subplot(3,1,1); plot(f3,20*log10(abs(CH./max(CH)))); grid; axis tight
  % xlabel('Frequency (MHz)')
  % ylabel('Power (dB)')
  % subplot(3,1,2); plot(uniform_win,'b'); grid;
  % axis([0 len 0 1.5])
  % subplot(3,1,3); plot(f3,20*log10(abs(chf_uni))); grid; axis tight
  % xlabel('Frequency (MHz)')
  % ylabel('Power (dB)')
  %
  %
  % % Apply triangular window to frequency-domain data and normalize
  % chf_tri = tri_win.*CH;
  % chf_tri = chf_tri./max(chf_tri);
  %
  % subplot(3,1,2); hold on; plot(tri_win,'r')
  % subplot(3,1,3); hold on; plot(f3,20*log10(abs(chf_tri)),'r'); grid; axis tight
  %
  %
  % % Apply hanning window to frequency-domain data and normalize
  % chf_hann = hann_win.*CH;
  % chf_hann = chf_hann./max(chf_hann);
  %
  % subplot(3,1,2); hold on; plot(hann_win,'g')
  % hold off
  % subplot(3,1,3); hold on; plot(f3,20*log10(abs(chf_hann)),'g'); grid; axis tight
  % hold off
  %
  % legend('Uniform','Triangle','Hanning')
  
  
  % Topic 4: Other Useful Functions
  
  % Interpolation
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
  
  fprintf('Section 9 Complete.\n');
  return
  
end



