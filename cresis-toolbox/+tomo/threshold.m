function twtt_sur_all = threshold(theta, img, Time ,slice)
% twtt_sur_all = threshold(theta, img, Time ,slice)
%
% Tract the ice-surface layer in the MUSIC slice.
%
% theta: theta contains 64 direction of arivals
% img: 3D imagery
% Time: two way travel time axis (twtt are points on the Time)
% slice: Integer from 1 to Nx (the currect working slice)
%
% twtt_sur_all: Two waytravel time (twtt) corresponding to each slice is
% stored in the respective column. It contains twtt of all the slices. 
%
% Author: Sravya Athinarapu
% Example:
% See slicetool_threshold.m



% index_nfloor limit choose properly
% for the datasets tested... for 039 Nfloor_limit = 0.63e-5;
% for 044 Nfloor_limit = 0.85e-5));
%similarly change it for other data sets

%  index_break limit choose properly --- this is used to separate the
%  surface return and the bottom return.


rline = slice;
num_pixels_col = length(Time);
pixel_search = floor(num_pixels_col./6);
 
% Request User Input

% figure,imagesc([],Time,squeeze(lp(img(:,:,rline))));
% prompt = 'enter the noise floor limit : \n ';
% Nfloor_limit = input(prompt);
% 
% prompt = 'separation between the surface return and the bottom return(for example 30 or 40 pixels separation) : \n ';
% separation= input(prompt);

 Nfloor_limit =  0.85e-5;
 separation = 40;

for idx = 1:length(theta)
    
    % IMAGE FILTER
    h1 = hann(20);
    h2 = hann(3);
    h = h1*h2.';
    img_filtered = filter2(h,squeeze(lp(img(:,:,rline))));
    rline_data = filter2(h,squeeze(lp(img(:,idx,rline))));
    index_nfloor = max (find(Time<= Nfloor_limit));
    time_thresholded = Time(index_nfloor:end);
    rline_data_thresholded = rline_data(index_nfloor:end);
    
    if idx == 1
        [req_pixel,id_req_pixel] = max(rline_data_thresholded);
        twtt_sur(idx) = time_thresholded(id_req_pixel);
    end
    
    if idx >= 2
        if idx <= 33
            data_DOA = rline_data_thresholded (1 :id_req_pixel);
            time_DOA = time_thresholded (1: id_req_pixel );
            corr = req_pixel.* data_DOA ;
            [req_pixel_both,id_req_pixel_both] = sort(corr,'descend');
            id_req_pixel = id_req_pixel_both(1);
            if idx <= 4                
                index_break =  id_req_pixel_both(min(find(abs(id_req_pixel_both - id_req_pixel)> separation)));
                id_req_pixel_both_sort = [index_break  id_req_pixel];
                [id_req_pixel_both_sorted,index] =  sort(id_req_pixel_both_sort,'ascend');
                req_pixel_both_sorted =req_pixel_both(index);
                twtt_s = time_DOA(id_req_pixel_both_sorted);
                twtt_sur(idx) = min(twtt_s);
            else                
                twtt_sur(idx) = time_DOA(id_req_pixel);
            end
        else   % for idx >34
            data_DOA = rline_data_thresholded;
            time_DOA = time_thresholded ;
            corr = req_pixel.* data_DOA ;
            [req_pixel_both,id_req_pixel_both] = sort(corr,'descend');
            id_req_pixel = id_req_pixel_both(1);
            index_break =  id_req_pixel_both(min(find(abs(id_req_pixel_both - id_req_pixel) > separation)));
            id_req_pixel_both_sort = [index_break  id_req_pixel];
            [id_req_pixel_both_sorted,index] =  sort(id_req_pixel_both_sort,'ascend');
            req_pixel_both_sorted =req_pixel_both(index);
            twtt_s = time_DOA(id_req_pixel_both_sorted);
            twtt_sur(idx) = min(twtt_s);
        end
    end
end

twtt_sur_all(:,rline) = twtt_sur';
return