function [ sources, references, num_slices] = setup_3Ddata(param, segments_and_frame)
% DATA_SETUP 
% Setup the data to run the tests
%
%     Input: 
%         segments_and_frame: a "hash map" that contains the segments and
%         the frames specific to each segment we'd like to test on
%
%     Output:
%          sources: array of filepaths(strings) of the frames that contains
%          the imagary
%          
%          references: array of references(strings) that contains the
%          ground truths with respect to its source
%
%          num_images: total number of 2D images from all the frames from the
%          sources
%
%          cpu_time: the estimated time that requires to compute the test
%          for each parameter test (test on all the frames)
% 
%          memory_requirement:  the estimated memory (Mbps) that requires
%          for each parameter test (test on all the frames)

    sources = {};
    references = {};

    for segment_name = keys(segments_and_frame)
      param.day_seg = segment_name{:};                                            % get the string
      
      for frame = segments_and_frame(segment_name{:})    
        source_fn = fullfile(ct_filename_out(param,param.out_type,''), ...
          sprintf('Data_%s_%03d.mat',param.day_seg,frame));                       % get the file path as a string from the image data
        ref_fn = fullfile(ct_filename_out(param,param.surfdata_source, ...
          'CSARP_surfData'),sprintf('Data_%s_%03d.mat',param.day_seg,frame));     % get the reference path as a string
        sources = [sources {source_fn}];
        references = [references {ref_fn}];
      end
    end

    num_slices = 0; % get the total number of slices (image) in the test
    for i = 1:length(sources)
      data=load(sources{i});  
      try
        images = data.Topography.img;
      catch ME
        error('Error in loading file: %s', ME.getReport);
      end  
      num_images_this_frame = size(images, 3);
      num_slices = num_slices + num_images_this_frame;
      clear images;
    end
    
end

