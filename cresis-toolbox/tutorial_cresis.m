%% script tutorial_cresis
%
% The purpose of this script is to demonstrate the use of several
% cresis functions and Matlab features that are relevant to processing
% cresis data.
%
% Author: John Paden 

lesson =1;

fprintf('\n\n########################## Lesson %d ##################\n\n', lesson);
%%
if lesson == 1
  % Cell arrays and get_filenames

  % To create an empty string
  filepath = '';
  
  % Create a file path string (Matlab data type char for character):
  if ispc
    % On Windows platform:
    filepath = 'X:\public\data\rds\2009_Greenland_TO\CSARP_standard\20090409_01\';
  else
    % On other platforms:
    filepath = '/cresis/snfs1/dataproducts/public/data/rds/2009_Greenland_TO/CSARP_standard/20090409_01/';
  end

  % Get the current working directory (pwd = print working directory)
  current_path = pwd;
  
  % Since filepath is a variable, we have to use the function form
  % for change directory (cd)
  cd(filepath);

  % Now list the contents of the directory for reference
  ls

  % Now go back to where we started
  cd(current_path);
  
  % Now use get_filenames to return all the properly matching files
  % in a cell array
  
  % Returns filenames matching "Data**.mat"
  filenames_all = get_filenames(filepath,'Data','','.mat');
  
  % Returns filenames satisfying "Data_20090409*.mat"
  %   In other words, the computer time on the radar system was:
  %    Year: 2009
  %    Month: 04
  %    Day: 09
  % UTC time: Universal Time Coordinate
  %   This is THE standard for time across the world. It is equivalent
  % to GMT or Greenwich Mean Time which is the time at the prime
  % meridian or zero deg longitude. GPS time is slightly different
  % than UTC time. Most people don't know this, but there are leap seconds
  % every few years to account for a small slow down in Earth's rotation
  % ... GPS time does NOT account
  % for these. However, the method for accounting for them is
  % fairly straight forward since the difference between GPS and UTC
  % is well documented and leap seconds (an extra second) in the UTC clock
  % so that it lags GPS time are always added at the end of the year.
  % For example, in 2009, GPS time is 15 seconds ahead of UTC.
  filenames_all = get_filenames(filepath,'Data_20090409','','.mat');

  % Get the first entry in filenames_all
  %   By indexing the cell array with {} I am extracting the contents
  %   of that cell array location.
  first_entry = filenames_all{1};
  
  % Create a 1 element cell array which has the first entry
  %   By indexing the cell array with () I am using the ordinary
  %   matrix indexing method. 
  first_entry_cell = filenames_all(1);
  
  % Useful short cuts:
  %  Instead of:
  last_entry = filenames_all{numel(filenames_all)};
  %  Use "end" to replace number of elements function numel():
  last_entry = filenames_all{end};

  % Create a cell array containing just the first three elements:
  entry1to3 = filenames_all(1:3);
  
  % Another trick which extracts the first three elements in the cell array
  [entry1 entry2 entry3] = filenames_all{1:3};
  
  % Search through each file name and pull out the frame indices in
  % the data filename. This is the last 3 digits of the filename.
  %   idx is short for index
  frms = [];
  for idx = 1:numel(filenames_all)
    % Useful Matlab command for getting just the file name without the
    % path to the file
    [path name ext] = fileparts(filenames_all{idx});
    % file_idx = extract last 3 characters from name and convert to double
    frms = [frms str2double(name(end-2:end))];
  end
  % Now create a new cell array with just the filenames associated with
  % frames 5 to frame 7
  %   Note: I have used () and not {}
  filenames_range = filenames_all(frms >= 5 & frms <= 7);
  
  % Just print out some examples by removing the semi-colon at the end of the line:
  filenames_range{1}
  filenames_range{end}
  
end
%%
if lesson == 2
  % Strings and printing
  
  % To put a ' in the middle of a string, simply type it twice:
  example_str = 'My example ''Here'' has two single quotes in it'
  
  % Print out the year month day in YYYYMMDD format
  year = 2009;
  month = 4;
  day = 9;
  fprintf('Normal:\n  %f/%f/%f\n', year, month, day);
  fprintf('Truncate to 1 decimal point:\n  %.1f/%.1f/%.1f\n', year, month, day);
  fprintf('Truncate to 0 decimal points:\n  %.0f/%.0f/%.0f\n', year, month, day);
  fprintf('Force minimum space:\n  %10f/%10f/%10f\n', year, month, day);
  fprintf('Force minimum space and 0 decimal points:\n  %4.0f/%2.0f/%2.0f\n', year, month, day);
  fprintf('Force minimum space with ''0'' and 0 decimal points:\n  %04.0f/%02.0f/%02.0f\n', year, month, day);
  fprintf('Take out slashes:\n  %04.0f%02.0f%02.0f\n', year, month, day);
  
  % Now create a YYYYMMDD string:
  ymd_str = sprintf('%04.0f%02.0f%02.0f', year, month, day);
  % Note create a pretend file name:
  filename = sprintf('data.%s080101.0000.raw', ymd_str)
  % Note create a pretend file path with Matlab's useful fullfile function
  % which puts in the correct slashes depending on the file system:
  filename = fullfile('tmp',ymd_str,'my_file')
end
 %% 
if lesson == 3
  %1). Write content into a specific file which is defined as the
  % 'filename' string variable
  filename = 'test_file_with_long_name.txt';
  
  % 'w': Create file for writing or overwrite if it already exists
  % note 1: for binary files need to watch out for 
  %   ieee little endian and ieee big endian formats
  % note 2: in Windows, text files have \r\n at the end and not
  %   just \n like Linux (\r = 13 = return, \n = 10 = newline)
  [fid,msg] = fopen(filename,'w');
  if fid == -1
    fprintf('Failed to open %s for writing\n', filename);
    fprintf(msg);
    return;
  else
    fprintf('Created file %s for writing\n', filename);
  end
  
  fprintf(fid, 'Here is the first value: 10\n');
  fprintf(fid, 'This is the second value: 11\n');
  fprintf(fid, 'Twelve is the third value: 12\n');

  fclose(fid);
  
  %2). Now read the file using textread
  %   This would not work if any one of the lines violates the format
  %   of 5 strings followed by 1 value. In that case we have to read
  %   each line in individually with textread and parse it.
  %   textread can use various delimiters and has lots of options...
  [S1 S2 S3 S4 S5 V1] = textread(filename,'%s %s %s %s %s %f');
  % Numbers, %f, are returned in a normal matrix
  index={'First','Second','Third'};
  for idx = 1:numel(V1)
    fprintf('%s value read: %f\n', index{idx},V1(idx));
  end
  % Strings, %s, are returned in a cell matrix
  for idx = 1:numel(V1)
    fprintf('%s string read: %s\n',index{idx}, S1{idx});
  end
end
%%
if lesson == 4
  % Matrix commands

  % Create a 3x3x2 matrix of normal gaussian variables
  A = randn([3 3 2])
  
  % Different ways to get A(1,3,2)
  % What this gets at is the matrix A is stored by marching down
  % rows/1st dimension, then columns/2nd dimension, then 3rd dim, 4th, etc
  % If you don't specify all the dimensions in a reference, all the
  % upper dimensions collapse down to the last entry
  row = 1;
  col = 3;
  third = 2;
  A(row,col,third)
  A(row,(third-1)*size(A,2) + col)
  A((third-1)*size(A,2)*size(A,1) + (col-1)*size(A,1) + row)
  
  % So to get all elements as one long vector (transposed with .')
  A1 = A(:).'
  % To get all row 1's and collapse the 2nd and 3rd dimensions:
  A2 = A(1,:)
  % To get all row 1's and keep structure to 2nd and 3rd dimensions
  A3 = A(1,:,:)
  % To just grab row 1 and third dimension 1
  A4 = A(1,:,1)
  % To just grab row 1 and third dimension 2
  A5 = A(1,:,2)
  % To just grab col 3 and third dimension 1
  A6 = A(:,3,1)
  
  % Get indices of all absolute values less than 0.1
  idxs_small = find(abs(A) < 0.3);
  A(idxs_small)
end
%%
if lesson == 5
  % Structures and checking
  
  clear param;
  
  fprintf('Exist %d\n', exist('param','var'));

  param = [];
  fprintf('Exist %d\n', exist('param','var'));
  fprintf('Struct %d\n', isstruct(param));
  
  % Create a structure array with no elements
  param = struct();
  fprintf('Struct %d\n', isstruct(param));
  
  param(1).field1 = 2;
  fprintf('field2 %d\n', isfield(param,'field2'));
  param(1).field2 = 3;
  fprintf('field2 %d\n', isfield(param,'field2'));
  
  param(1).field1 = 4
  % This is a special case which only works when param is a 1 by 1
  % struct array
  param.field1 = 4

  param(2).field1 = 5
  param(1).field1
  param(2).field1
  try
    param.field1 = 6;
  catch
    fprintf('Did not work since param is a 2x1 struct array now\n');
    fprintf('Entering debug mode (i.e. hard breakpoint)\n');
    fprintf('Commands: dbcont, dbstep, dbstep in, dbstep out, dbquit, dbstack\n');
    fprintf('Type dbcont to continue or dbquit to quit\n');
    fprintf('  or use debug commands in editor\n');
    keyboard
  end
  
end
%%
if lesson == 6
  % Create a 3x3 matrix of complex normal gaussian variables
  % where variance is 1
  noise = 1/sqrt(2) * (randn([3 3]) + j*randn([3 3]))

  % Create an Nt by Nx matrix of complex normal gaussian variables
  % where variance is 1
  % Nt = number of fast time samples (first dimension)
  Nt = 600;
  % Nx = number of spatial samples (rangelines, records, slow time, second
  % dimension)
  Nx = 1000;
  noise = 1/sqrt(2) * (randn([Nt Nx]) + j*randn([Nt Nx]));

  % Array of all dimensions of noise
  size(noise)
  % Scalar of first dimension of noise
  size(noise,1)
  % Scalar of second dimension of noise
  size(noise,2)
  
  % Check that the mean variance is nearly 1
  mean_noise_pow = mean(abs(noise(:)).^2)
  
  % Filter the first range line of the data with 2nd order butterworth IIR
  % filter and a cutoff of 0.02 * fs/2 where fs = sampling frequency
  %   B is the feed forward (FIR) coefficients
  %   A is the feed back (IIR) coefficients, A = 1 for FIR
  %   IIR = infinite impulse response
  %   FIR = finite impulse response
  line = 1;
  [B,A] = butter(2,0.02);
  noise_row1_filt = filter(B,A,noise(:,line));
  
  % Plot original and filtered in figure 1:
  figure(1); 
  plot(10*log10(noise(:,line) .* conj(noise(:,line))));
  hold on;
%   handle = plot(10*log10(noise_row1_filt .* conj(noise_row1_filt)),'r');
  plot(10*log10(noise_row1_filt .* conj(noise_row1_filt)),'r','LineWidth',2);
  hold off;
  
  % Make the filtered line thicker
%   set(handle,'LineWidth',2);
  
  grid on;
  
  xlabel('Range bins');
  ylabel('Relative power (dB)');
  title(sprintf('Original and filtered bins for line %d', line));
  
  % Use tex interpreter to create interesting axes labels
  %   tex interpreter is very powerful and can represent matrices, greek
  %   letters, subscripts, etc
  %   Sometimes lack of font support breaks it though (e.g. X-windows)
  ylabel('Relative power (10*log10(|noise|^2))');

  % Create a time axis
  % fs = sampling frequency
  fs = 120e6;
  % dt = sample time spacing
  dt = 1/fs;
  % time = fast time axis
  time = (0:dt:(Nt-1)*dt).';
  
  figure(2); 
  plot(time*1e6, 10*log10(noise(:,line) .* conj(noise(:,line))));
  hold on;
  handle = plot(time*1e6, 10*log10(noise_row1_filt .* conj(noise_row1_filt)),'r');
  hold off;
  grid on;
  
  % Should be microseconds (us), but might show up as (ms) if tex breaks
  xlabel('Time ({\mu}s)');
  ylabel('Relative power (dB)');
  
  % Remove delay and plot
  line = 1;
  [B,A] = butter(2,0.02);
  % Zero delay filtering (does some fancy tricks to create zero delay
  % without looking bad on the edges... usually does a good job with
  % edge conditions, but not always)
  noise_row1_filt = filtfilt(B,A,noise(:,line));
  
  figure(3); 
  plot(time*1e6, 10*log10(noise(:,line) .* conj(noise(:,line))));
  hold on;
  plot(time*1e6, 10*log10(noise_row1_filt .* conj(noise_row1_filt)),'r');
  hold off;
  grid on;
  xlabel('Time ({\mu}s)');
  ylabel('Relative power (dB)');
  
  % Show image
  figure(4); 
  h=imagesc(20*log10(abs(noise)));
 
  colormap(jet(256));
  xlabel('Slow time (lines)');
  ylabel('Fast time (bins)');
  handle = colorbar;
  set(get(handle,'YLabel'), 'String', ' Relative  Power  (dB)');


  %Try gray map
  % colormap(gray(256));
  
  % Close figure
  close(1);
end

% imagesc() function basic understanding
%%
if lesson==7
  % suppose we load data from the
  % \cresis\scratch2\mdce\accum\2011_Greenland_P3
  
  filepath ='Z:\mdce\accum\2011_Greenland_P3\CSARP_qlook\20110329_04\Data_20110329_04_009';
  %load data from this file, it will return a sturcture 'data';
  data=load(filepath),
  %get the time range
  Time=data.Time;
  Data=data.Data;
  figure(1)
 
  imagesc([],Time,lp(Data));
  %some other useful syntax for "imagesc" function
  %1. display the image in a certain range of x- and y-axis
  %   imagesc(x,y,C), where 'x' and 'y' is a vector which is used to specify
  %   the image bound on x- and y-axis. 'C' is a M-by-M matrix that each element
  %   represent a pixel in the colormap.
  %2. specify the display range either on the x-axis or on the y-axis.
  %   imagesc([],y,C) or imagesc(x,[],C), the blank bracket '[]' indicates a
  %   default display bound on x- or y-axis.
  
  % Try the following selected bound colormap display;
  %    imagesc([0 500],[0 300],lp(Data));
  %    imagesc([],[0 300],lp(Data));
  %    imagesc([0 500],[],lp(Data));

  
 
   %colormap(gray);
  
end
  
  
 %% 
% Plot handle manipulation 
 if lesson == 8
   x = -pi:.1:pi;
   y = sin(x);
   z = cos(x);
   %p is the handle of this plot
   p = plot(x,y);
   
   
   get(p);
   %viwe the various properties of the plot
   
   %using get(handle,...) and set(handle,...) functions to view the plot
   %properties and set a certain plot property, respectively.
   %gca, get current axes properties. It will return a handle of the current axis
   %gcf, get current figure properties. It will return a handle of the current figure.
   %get(0), returns the current values of all user-settable properties.
   %for example, if we want to display the range of the x-axis from -pi to pi with space
   %-pi/2
   %the 'XTick'characteristic inside the axes handle can be set as
   set(gca,'XTick',-pi:pi/2:pi);
   
   %In addition, we also want to label the x-axis as the the setup we defined
   %in the last step, then we can use 'XTickLabel' inside the handle 'p'
   set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
   %Finally, we can use xlabel() and ylabel() function to set the display name
   %of the x- and y- axis
   xlabel('-\pi \leq \Theta \leq \pi')
   ylabel('sin(\Theta)')
   title('Plot of sin(\Theta)')
   %*more details about various characteristics of axes is availabe on the website
   %http://www.mathworks.com/help/techdoc/ref/axes_props.html
   
   
   %there is another useful feature which is 'Visible' inside of the handle
   %property. 'Visible' has two status which are 'on' and 'off'. when the
   %status is 'off' the plot will not display on the figure, but the handle of
   %this plot is still availabe
   %for example
   set(p,'Visible','off')
   pause; %press anykey in command window will continue the execution
   set(p,'Visible','on')
   
   %the advantage of utilizing 'Visible' is that it can save more time
   %to display the desired plot than using two plot() functions
   %the comparsion is shwon below: (if the finall plot is cos(x) )
   %1. using two plot() finctions
   figure(2)
   start1=tic;
   plot(x,y);
   plot(x,z);
   elapse1=toc(start1)
   
   figure(3)
   start2=tic;
   set(p,'Visible','off');
   plot(x,z);
   elapse2=toc(start2)
 
   
   %Comment: The Figure(2) and Figure(3) have the same display result,
   %however, you could find that the elapse1>elapse2. It's indicate that by
   %applying the 'Visible' property, the execution time is shorter than directly
   %update the plot with two plot() functions. Especially, when the plot is
   %dealing with a large data, the significance of 'Visible' will be overwhelmed
   
 end
 
return;
%%