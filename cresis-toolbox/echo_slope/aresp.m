function aresp(dir_block, file_block, decim_x, decim_z, z_smooth, num_overlap, parallel_check, area_max, ellipse_min, angle_max, dist_max, num_obj_max)
% ARESP Calculate layer slopes using image-processing techniques.
%   
%   ARESP uses image-processing techniques to transfrom a radargram into a
%   binary image, filter it, and then extract layer-slope information from
%   that image. In this implementation, radargram "blocks" outputted by
%   RADBLOCKPROC are loaded, analyzed and then saved with an additional
%   field, block.slope_aresp, which is the inferred layer slope at each
%   point in the radargram. This layer slope can be used by PICKGUI to
%   predict the radiostratigraphy and then flatten the radargram.
%   
%   ARESP(DIR_BLOCK,FILE_BLOCK,DECIM_X,DECIM_Z,Z_SMOOTH,NUM_OVERLAP,PARALLEL_CHECK,AREA_MAX,ELLIPSE_MIN,ANGLE_MAX,DIST_MAX,NUM_OBJ_MAX)
%   performs the operations described above. DIR_BLOCK is the directory
%   containing the radargram "blocks". FILE_BLOCK is the filename pattern
%   (accepts wildcards) of blocks to be processed. DECIM_X is a two-element
%   vector, where each element is the number of indices to smooth (L_x) and
%   decimate (2) the binary radargram by in the along-track direction.
%   DECIM_Z is the number of indices to decimate the binary radargram by in
%   the vertical direction. Z_SMOOTH is a two-element vector, where each
%   element is a number of indices by which two binary versions of the
%   original radargram are smoothed vertically (L_z1 and L_z2). NUM_OVERLAP
%   is a two-element vector, where each element is the number of overlaps
%   between vertical slices of the two versions of the binary radargram to
%   be examined by ARESP_OBJ. PARALLEL_CHECK is a logical scalar that is
%   true (1) if a license for the Parallel Computing Toolbox is available,
%   in which case certain operations will be parallelized. AREA_MAX is the
%   maximum area of a valid layer object (A_2). ELLIPSE_MIN is the minimum
%   ellipticity of a valid layer object (r). ANGLE_MAX is the maximum valid
%   slope of a layer object. DIST_MAX is the maximum distance in index
%   space to consider a fit between the decimated binary radargram and the
%   centroid of a layer object. NUM_OBJ_MAX is the maximum number of layer
%   objects used to calculate the median object slope at any given point in
%   the decimated binary radargram.
%   
%   ARESP was originally developed and described by:
%   
%   Sime, L.C., R.C.A. Hindmarsh and H.F.J. Corr (2011), Automated
%   processing to derive dip angles of englacial radar reflectors in ice
%   sheets, J. Glaciol., 57(202), 260-266.
%   
%   The notation used above in parentheses is from this paper.
%   
%   ARESP is based on functions made available at:
%   
%   http://researchpages.net/QDES/people/louise-sime/resources/
%   
%   See also ARESP_OBJ.
%   
% Louise Sime (BAS), Joe MacGregor (UTIG)
% Last updated: 06/27/14

% clear
% % dir_block                   = '/Volumes/icebridge/data/2006_to/block/20060601_01/';
% dir_block                   = '~/Desktop/';
% file_block                  = '20130402_01_block10_43';
% decim_x                     = [10 10];
% decim_z                     = 5;
% z_smooth                    = [6 36];
% num_overlap                 = [2 2];
% parallel_check              = false;
% area_max                    = 3; % maximum ratio of object area to ind_sep*z_smooth
% ellipse_min                 = 2; % minimum object ellipticity (MajorAxisLength/MinorAxisLength)
% angle_max                   = 30; % maximum object orientation
% dist_max                    = 25; % maximum number of indices
% num_obj_max                 = 300; % maximum number of objects
% nargin                      = 12;

if ~license('test', 'image_toolbox')
    error('aresp:image', 'Image Processing Toolbox license is not available.')
end
if (nargin ~= 12)
    error('aresp:nargin', 'Not enough input arguments (need 12).')
end
if ~ischar(dir_block)
    error('aresp:dirblockstr', 'Block directory (DIR_BLOCK) is not a string.')
end
if (~isempty(dir_block) && ~exist(dir_block, 'dir'))
    error('aresp:nodirblock', 'Block directory (DIR_BLOCK) does not exist.')
end
if ~ischar(file_block)
    error('aresp:fileblockstr', 'Block filename (FILE_BLOCK) is not a string.')
end
if (~isnumeric(decim_x) || ~isvector(decim_x) || (length(decim_x) ~= 2))
    error('aresp:decimxtype', 'Horizontal decimation (DECIM_X) is not a two-element numeric vector.')
end
if any(mod(decim_x, 1))
    decim_x                 = round(decim_x);
    warning('aresp:rounddecimx', ['Horizontal decimation (DECIM_X) rounded to [' num2str(decim_x) '].'])
end
if (~isnumeric(decim_z) || ~isscalar(decim_z))
    error('aresp:decimztype', 'Vertical decimation (DECIM_Z) is not a numeric scalar.')
end
if mod(decim_z, 1)
    decim_z                 = round(decim_z);
    warning('aresp:rounddecimz', ['Vertical decimation (DECIM_Z) rounded to ' num2str(decim_z) '.'])
end
if (~isnumeric(z_smooth) || ~isvector(z_smooth) || (length(z_smooth) ~= 2))
    error('aresp:zsmoothtype', 'Vertical filter size (Z_SMOOTH) is not a two-element numeric vector.')
end
if any(mod(z_smooth, 1))
    z_smooth                = round(z_smooth);
    warning('aresp:roundzsmooth', ['Vertical filter size (Z_SMOOTH) rounded to [' num2str(z_smooth) '].'])
end
if (~isnumeric(num_overlap) || ~isvector(num_overlap) || (length(num_overlap) ~= 2))
    error('aresp:numoverlaptype', 'Vertical filter size (NUM_OVERLAP) is not a two-element numeric vector.')
end
if any(mod(num_overlap, 1))
    num_overlap             = round(num_overlap);
    warning('aresp:roundnumoverlap', ['Vertical filter size (NUM_OVERLAP) rounded to [' num2str(num_overlap) '].'])
end
if (~islogical(parallel_check) || ~isscalar(parallel_check))
    error('aresp_obj:parallelcheck', 'True/false check for Parallel Computing Toolbox license (PARALLEL_CHECK) is not a logical scalar.')
end
if (~isnumeric(dist_max) || ~isscalar(dist_max))
    error('aresp:distmaxtype', 'Maximum distance (DIST_MAX) is not a numeric scalar.')
end
if (~isnumeric(num_obj_max) || ~isscalar(num_obj_max))
    error('aresp:numobjmaxtype', 'Maximum number of objects (NUM_OBJ_MAX) is not a numeric scalar.')
end
if mod(num_obj_max, 1)
    num_obj_max             = round(num_obj_max);
    warning('aresp:roundnumobjmax', ['Maximum number of objects (NUM_OBJ_MAX) rounded to ' num2str(num_obj_max) '.'])
end
if (~isnumeric(area_max) || ~isscalar(area_max))
    error('aresp:areamaxtype', 'Maximum object area (AREA_MAX) is not a numeric scalar.')
end
if (~isnumeric(ellipse_min) || ~isscalar(ellipse_min))
    error('aresp:ellipsemintype', 'Minimum object ellipticity (ELLIPSE_MIN) is not a numeric scalar.')
end
if (~isnumeric(angle_max) || ~isscalar(angle_max))
    error('aresp:anglemaxtype', 'Maximum object angle (ANGLE_MAX) is not a numeric scalar.')
end

ind_sep                     = [(z_smooth(1) * 3) (z_smooth(1) * 6)]; % S_1, S_2

% filenames in dir_block
name_file                   = dir([dir_block file_block '.mat']);
name_file                   = {name_file.name}';
num_file                    = length(name_file);

disp(['Starting ARESP for ' dir_block '...'])

% loop through each file in dir_block
for ii = 1:num_file
    
    % load current file
    disp([name_file{ii}(1:(end - 4)) ' (' num2str(ii) ' / ' num2str(num_file) ')...'])
    block                   = load([dir_block name_file{ii}], 'block');
    block                   = block.block;
    
    % remove old version of field if present
    if isfield(block, 'slope_aresp')
        block               = rmfield(block, 'slope_aresp');
    end
    
    % trim just below bed if possible
    if isfield(block, 'twtt_bed')
        try
            num_sample_trim = interp1(block.twtt, 1:block.num_sample, nanmax(block.twtt_bed), 'nearest', 'extrap');
        catch %#ok<CTCH>
            num_sample_trim = block.num_sample;
        end
    else
        num_sample_trim     = block.num_sample;
    end
    if isnan(num_sample_trim)
        num_sample_trim     = block.num_sample;
    end
    amp_mean                = block.amp(1:num_sample_trim, :);
    
    % horizontal smoothing
    tmp1                    = floor(decim_x(1) / 2);
    for jj = (1 + ceil(decim_x(1) / 2)):(block.num_trace - ceil(decim_x(1) / 2))
        amp_mean(:, jj)     = nanmean(block.amp(1:num_sample_trim, (jj - tmp1):(jj + tmp1)), 2); % P_x
    end
       
    % vertical smoothing
    try
        amp_filt_hi         = colfilt(amp_mean, [z_smooth(1) 1], 'sliding', @nanmean); % hi-res z-filtering, P_xz1
        amp_filt_lo         = colfilt(amp_mean, [z_smooth(2) 1], 'sliding', @nanmean); % lo-res z-filtering, P_xz2
    catch %#ok<CTCH>
        disp('Vertical filtering failed on this block...')
        continue
    end
    
    % use averaged arrays to generate binary arrays
    amp_bin_hi              = (amp_mean - amp_filt_lo) < 0; % hi-res binary radargram, B_1
    amp_bin_lo              = (amp_filt_hi - amp_filt_lo) < 0; % lo-res binary radargram, B_2
    
    % object properties
    bin_hi_obj              = aresp_obj(amp_bin_hi, ind_sep(1), num_overlap(1), parallel_check); % hi-res objects
    bin_lo_obj              = aresp_obj(amp_bin_lo, ind_sep(2), num_overlap(2), parallel_check); % lo-res objects
    
    % exclude invalid object data
    bin_hi_obj              = bin_hi_obj(((bin_hi_obj(:, 3) <= (ind_sep(1) * z_smooth(1) * area_max)) & (bin_hi_obj(:, 3) >= ind_sep(1)) & ((bin_hi_obj(:, 4) ./ bin_hi_obj(:, 5)) >= ellipse_min) & (abs(bin_hi_obj(:, 6)) <= angle_max)), [1 2 6]);
    bin_lo_obj              = bin_lo_obj(((bin_lo_obj(:, 3) <= (ind_sep(2) * z_smooth(2) * area_max)) & (bin_lo_obj(:, 3) >= ind_sep(2)) & ((bin_lo_obj(:, 4) ./ bin_lo_obj(:, 5)) >= ellipse_min) & (abs(bin_lo_obj(:, 6)) <= angle_max)), [1 2 6]);
    
    % all valid objects
    obj_good                = [bin_hi_obj; bin_lo_obj];
    
    % horizontal decimation
    ind_decim_x             = [1 (1 + ceil(decim_x(2) / 2)):decim_x(2):(block.num_trace - ceil(decim_x(2) / 2)) block.num_trace];
    num_decim_x             = length(ind_decim_x);
    
    % vertical decimation
    ind_decim_z             = ((1 + ceil(decim_z / 2)):decim_z:(num_sample_trim - ceil(decim_z / 2)))';
    num_decim_z             = length(ind_decim_z);
    
    % calculate median object orientation/angle
    slope_med               = zeros(num_decim_z, num_decim_x);
    if parallel_check
        [tmp1, tmp2, tmp3]  = deal(obj_good(:, 1), obj_good(:, 2), obj_good(:, 3)); % slice for parfor
        parfor jj = 1:num_decim_z
            for kk = 1:num_decim_x
                dist        = sqrt(((ind_decim_x(kk) - tmp1) .^ 2) + (((ind_decim_z(jj) - tmp2) .^ 2))); %#ok<PFBNS> % distance
                ind_close   = find(dist < dist_max);
                if (length(ind_close) > num_obj_max) % if more than num_obj_min objects are nearby, use only the closest ones
                    [~, ind_close_sort] ...
                            = sort(dist(ind_close)); % closest objects
                    slope_med(jj, kk) ...
                            = median(tmp3(ind_close(ind_close_sort(1:num_obj_max)))); %#ok<PFBNS>
                elseif ~isempty(ind_close) % otherwise use all
                    slope_med(jj, kk) ...
                            = median(tmp3(ind_close));
                end
            end
        end
    else
        for jj = 1:num_decim_z
            for kk = 1:num_decim_x
                dist        = sqrt(((ind_decim_x(kk) - obj_good(:, 1)) .^ 2) + ((ind_decim_z(jj) - obj_good(:, 2)) .^ 2)); % distance
                ind_close   = find(dist < dist_max);
                if (length(ind_close) > num_obj_max)
                    [~, ind_close_sort] ...
                            = sort(dist(ind_close));
                    slope_med(jj, kk) ...
                            = median(obj_good(ind_close(ind_close_sort(1:num_obj_max)), 3));
                elseif ~isempty(ind_close)
                    slope_med(jj, kk) ...
                            = median(obj_good(ind_close, 3));
                end
            end
        end
    end
    
    [tmp1, tmp2]            = meshgrid(ind_decim_x, ind_decim_z);
    [tmp3, tmp4]            = meshgrid(1:block.num_trace, 1:num_sample_trim);
    block.slope_aresp       = interp2(tmp1, tmp2, tand(slope_med), tmp3, tmp4, '*spline', 0);
    if (num_sample_trim < block.num_sample)
        block.slope_aresp((num_sample_trim + 1):block.num_sample, :) ...
                            = NaN;
    end
    block                   = orderfields(block); %#ok<NASGU>
    save([dir_block name_file{ii}], '-v7.3', 'block')
    
end

% figure
% subplot(211)
% imagesc(amp_bin_hi)
% colormap(jet)
% colorbar('fontsize', 20)
% title(['Binary ' name_file{ii}(1:(end - 4))], 'interpreter', 'none')
% subplot(212)
% imagesc(ind_decim_x, ind_decim_z, slope_med, [-angle_max angle_max])
% colorbar('fontsize', 20)
% title('ARESP angle')