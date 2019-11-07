% get_google_map returns a Google Static Map depending on the map_zone
% selected and centered at the default lat lon set for the map_zone
function A = get_google_map(obj)

  %% Setting default lat lon depending on the map_zone
  if strcmpi('arctic', obj.map_pref.settings.map_zone)
    c_lat = 73.82177;
    c_lon = -40.333279;
    zoom = 3;
  else
    % Antarctica
    c_lat = -73.381418;
    c_lon = -67.900671;
    zoom = 3;
  end
  
  % Creating a Picker Google Map object
  googleObj = imb.picker_google_map;
  
  % Requesting a Google Staic Map from the server
  googleObj = imb.request_google_map(googleObj, c_lat, c_lon, zoom);
  
  % Updating the map object
  obj.googleObj = googleObj;
  
  % Returning the Map
  A = obj.googleObj.A;
end