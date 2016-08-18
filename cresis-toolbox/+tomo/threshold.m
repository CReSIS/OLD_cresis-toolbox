function twtt_sur_all = threshold(mdata,slice)
% twtt_sur_all = threshold(mdata,slice)
%
% Tract the ice-surface layer in the MUSIC slice.
%
% mdata:
% slice:
%
% twtt_sur_all:
%
% Author: Sravya Athinarapu

% index_nfloor limit choose properly
% for the datasets tested... for 039 index_nfloor = max (find(mdata.Time<= 0.63e-5));
% for 044 index_nfloor = max (find(mdata.Time<= 0.85e-5)); 
%similarly change it for other data sets
%  index_break limit choose properly

DOA_used = mdata.param_combine.array_param.theta*(180/pi);
rline = slice;
num_pixels_col = length(mdata.Time);
pixel_search = floor(num_pixels_col./6);

for DOA = 1:64
    
    % IMAGE FILTER    
    h1 = hann(20);
    h2 = hann(3);
    h = h1*h2.';
    img_filtered = filter2(h,squeeze(lp(mdata.Topography.img(:,:,rline))));    
    rline_data =   filter2(h,squeeze(lp(mdata.Topography.img(:,DOA,rline))));
    
    % change the limit if required  (depends for a dataset)
    index_nfloor = max (find(mdata.Time<= 0.63e-5));
    %index_nfloor = max (find(mdata.Time<= 0.85e-5));
    time_thresholded = mdata.Time(index_nfloor:end);
    rline_data_thresholded = rline_data(index_nfloor:end);
    
    if DOA == 1
        [req_pixel,id_req_pixel] = max(rline_data_thresholded);
        twtt_sur(DOA) = time_thresholded(id_req_pixel);
    end
    
    if DOA >=2 
        if DOA <=33            
            data_DOA = rline_data_thresholded (1 :id_req_pixel);
            time_DOA = time_thresholded (1: id_req_pixel );
            corr = req_pixel.* data_DOA ; 
            [req_pixel_both,id_req_pixel_both] = sort(corr,'descend');
                id_req_pixel = id_req_pixel_both(1);
            if DOA<= 4  
               
                index_break =  id_req_pixel_both(min(find(abs(id_req_pixel_both - id_req_pixel)> 40)));                
                id_req_pixel_both_sort = [index_break  id_req_pixel];                
                [id_req_pixel_both_sorted,index] =  sort(id_req_pixel_both_sort,'ascend');                
                req_pixel_both_sorted =req_pixel_both(index);                
                twtt_s = time_DOA(id_req_pixel_both_sorted);                
                twtt_sur(DOA) = min(twtt_s);                
            else 
                                
                twtt_sur(DOA) = time_DOA(id_req_pixel);                
            end                        
        else   % for DOA >34
            data_DOA = rline_data_thresholded;
            time_DOA = time_thresholded ;         
          corr = req_pixel.* data_DOA ;            
            [req_pixel_both,id_req_pixel_both] = sort(corr,'descend');
            id_req_pixel = id_req_pixel_both(1);            
            index_break =  id_req_pixel_both(min(find(abs(id_req_pixel_both - id_req_pixel)> 40)));            
            id_req_pixel_both_sort = [index_break  id_req_pixel];            
            [id_req_pixel_both_sorted,index] =  sort(id_req_pixel_both_sort,'ascend');            
            req_pixel_both_sorted =req_pixel_both(index);            
            twtt_s = time_DOA(id_req_pixel_both_sorted);            
            twtt_sur(DOA) = min(twtt_s);            
        end
    end
end

twtt_sur_all(:,rline) = twtt_sur';
return