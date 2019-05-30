classdef picker_google_map
  properties
    key = 'AIzaSyCNexiP6WcIda8ZEa2MnwznWrGotDoLu0w'
    tile_size = 256
    w = 1280
    h = 1280
    
    zoom % Zoom level
    A % Map image
    
    % Center coordinates
    c_lat 
    c_lon
    c_px_x
    c_px_y
    c_wc_x
    c_wc_y
    
    % Pixel coordinates of vertices
    top_left_px_x
    top_left_px_y
    bottom_left_px_x
    bottom_left_px_y
    top_right_px_x
    top_right_px_y
    bottom_right_px_x
    bottom_right_px_y
    
    % World Coordinates of vertices
    top_left_wc_x
    top_left_wc_y
    bottom_left_wc_x
    bottom_left_wc_y
    top_right_wc_x
    top_right_wc_y
    bottom_right_wc_x
    bottom_right_wc_y
    
    % Toggle values to see if user has requested to zoom in/out or not
    zoom_in_out
    
    % Stores the arrow key pressed by user to pan
    pan

    hfig
    hImage
  end
  
  %% README
  %
  % px: pixel coordinates
  % wc: world coordinates
  % c_*: center coordinates
  % pan: 1 - up arrow
  %      2 - right arrow
  %      3 - down arrow
  %      4 - left arrow
  
end