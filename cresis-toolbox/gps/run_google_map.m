global gGoogle;
try
  delete(gGoogle);
end

gGoogle = google_map();

% [A,x_axis,y_axis] = gGoogle.request_google_map(10,23,10,23);

% [A,x_axis,y_axis] = gGoogle.request_google_map(125,131,90,95);
[A,x_axis,y_axis] = gGoogle.greenland();
% [A,x_axis,y_axis] = gGoogle.antarctica();

% [A,x_axis,y_axis] = gGoogle.request_google_map(0,256,0,256);

figure(1); clf;
imagesc(x_axis,y_axis,A)
