%% first create the folders 
addpath('helper_functions');
create_folders; 
% first run these two lines to create folders under /script that will contain your test results
% this will only need to run once

% *IMPORTANT* PLEASE RUN EVERYTHING UNDER THIS FOLDER TO MAKE SURE THINGS ARE SAVED IN
% THE RIGHT PLACE.

%% PARAMETER TUNING 

% These files are used for tuning the parameters for Viterbi and TRWS 
% algorithms on both 3D and 2D data so that we can obtain better ice-bottom 
% detection to ease the burden of manual data labeling.

% tuning_method: Either Grid Search or Random Search
% algorithm: Viterbi, TRWS
% dataType: 2D or 3D 

% Ideally, the only file that needs to be edited are the
% user_settings_<tuning_method>_<algorithm>_<dataType>.m 

% After the user_settings_<algorithm>_<method>_<dataType>.m is set, go to the
% run_<tuning_method>_<algorithm>_<dataType>.m and run the script

% Please run the script run_<tuning_method>_<algorithm>_<dataType>.m under
% its folder so the result can be saved into its designated folder

%% How to tune the parameters

% Each file's name shows the tuning method, algorithm(TRWS or Viterbi), and the type of data
% (3D or 2D) to tested on.

% 1) Check the user_settings_<tuning_method>_<algorithm>_<dataType>.m in
% each folder.

% For 3D data, we use a Hash map to document the segments and the frames. The
% 'keys' for the Hash map are the segments, and the 'values' for each key
% are the frames. Edit those for the frames you'd like to run the test on.

% For 2D data, it simply refers to the filename of the param spreadsheet.

% For parameters' setup in Grid Search, you want to manually insert the points for each of the
% parameter for you'd like to test. The program will search the Cartesian product 
% of the user-defined subset of the range of each parameter. 

% For parameters' setup Random Search, you want to define the upper bound and the lower
% bound for each parameter. The program will generate random vectors based
% on this.

% 2) Go to run_<tuning_method>_<algorithm>_<dataType>.m and press 'F5' on the keyboard

% 3) The result will be automatically saved under the /scripts folder where
% it also contains the cresis-toolbox. 


%% Flow of the program

% 1) User setup the settings in the user_settings_<tuning_method>_<algorithm>_<dataType>.m
% 2) Run the script to start grid search or random search 
% 3) When the search is complete, the best_parameter (ranked via mean
% error) will be set as the initial point in fmincon. The subroutine
% fmincon will search the region near the initial point to find the optimal
% point in the parameter space using the same dataset we've used as in the
% search.
% 4) When the fmincon is done, the result will be saved.

% a) Repeat the above a couple more times using different dataset (different frames in different segment in 3D, for example). 
% b) We can then run the whole season by using the results obtained from these by running run_whole_season_<algorithm>_<datatype>.m 
% to see which combination of parameter is optimal.


%% MORE

% To know more about random search, please refer to:
% http://www.jmlr.org/papers/volume13/bergstra12a/bergstra12a.pdf
% Or just google "random search bengio"

%% NOTES

% *****IF YOU RUN 2D tests ON ANTARTICA**********, please make sure
% cluster_kernel_Viterbi_2D.m line 174 to change (if 0) to (if 1). 
% If RUN on GREENLAND then change line 174 from (if 1) to (if 0). 

% For some reason, the data in 2014032505 doesn't work; do not load it into
% Viterbi 3D tests

% The entrire season for 3D file is 2014 Greenland P3
% The entrire season for 2D file is 2009 Antartica DC8

% You can setup the number of iteration of fmincon. Just go to
% result_find_optimal_param_<dataType>.m line 18 to adjust it.


%% SHANE's contact

% guixien@gmail.com
% 917-640-7301







