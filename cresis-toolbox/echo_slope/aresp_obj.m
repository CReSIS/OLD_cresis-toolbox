function obj                = aresp_obj(data_bin, ind_sep, num_overlap, parallel_check)
% ARESP_OBJ Measure key object properties of binary radargrams.
% 
%   OBJ = ARESP_OBJ(DATA_BIN,IND_SEP,NUM_OVERLAP,PARALLEL_CHECK) separates
%   binary radargram DATA_BIN (B_1 or B_2) into vertical slices using a
%   specified index separation IND_SEP (S_1 or S_2) and number of overlaps
%   NUM_OVERLAP and returns OBJ, which contains the key object properties
%   of all slices. The same properties are also measured on the logical
%   inverse of DATA_BIN (~DATA_BIN). OBJ is an M x 6 array, where M is the
%   number of measured objects and 6 is the number of necessary object
%   properties plus one, because the Centroid property returns two values.
%   
%   If a license for the Parallel Computing Toolbox is available
%   (PARALLEL_CHECK = true), then ARESP_OBJ will be parallelized.
%   
%   Standard object properties for ARESP are: Centroid, Area,
%   MajorAxisLength, MinorAxisLength, Orientation. See REGIONPROPS
%   documentation for complete object property list.
%   
%   ARESP_OBJ was originally developed and described by:
%   
%   Sime, L.C., R.C.A. Hindmarsh and H.F.J. Corr (2011), Automated
%   processing to derive dip angles of englacial radar reflectors in ice
%   sheets, J. Glaciol., 57(202), 260-266.
%   
%   The notation used above in parentheses is from this paper.
%   
%   ARESP_OBJ is based on functions made available at:
%   
%   http://researchpages.net/QDES/people/louise-sime/resources/
%   
%   See also ARESP.
%   
% Louise Sime (BAS), Joe MacGregor (UTIG)
% Last updated: 06/27/14

if (nargin ~= 4)
    error('aresp_obj:nargin', 'Not enough input arguments (need 4).')
end
if ~islogical(data_bin)
    error('aresp_obj:databintype', 'Binary radargram array (DATA_BIN) is not a logical array.')
end
if (~isnumeric(ind_sep) || ~isscalar(ind_sep) || mod(ind_sep(1), 1))
    error('aresp_obj:indsep', 'Index separation (IND_SEP) is not a round numeric scalar.')
end
if (~isnumeric(num_overlap) || ~isscalar(num_overlap) || mod(num_overlap(1), 1))
    error('aresp_obj:indsep', 'Number of overlaps (NUM_OVERLAP) is not a round numeric scalar.')
end
if (~islogical(parallel_check) || ~isscalar(parallel_check))
    error('aresp_obj:parallelcheck', 'True/false check for Parallel Computing Toolbox license (PARALLEL_CHECK) is not a logical scalar.')
end

ind_step                    = ceil(ind_sep ./ num_overlap);
obj                         = [];

for ii = 1:2
    if (ii == 2) % also do logical inverse
        data_bin            = ~data_bin;
    end
    if parallel_check
        obj_tmp             = cell(num_overlap, 1);
        parfor jj = 1:num_overlap
            % slice binary radargram
            data_slice      = data_bin;
            data_slice(:, (jj * ind_step):ind_sep:end) ...
                            = false;
            % measure object properties of cleaned, sliced binary radargram
            curr_obj        = regionprops(labelmatrix(bwconncomp(bwmorph(data_slice, 'clean'), 4)), 'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
            % concatenate properties
            obj_tmp{jj}     = [cat(1, curr_obj.Centroid) [curr_obj(:).Area]' [curr_obj(:).MajorAxisLength]' [curr_obj(:).MinorAxisLength]' [curr_obj(:).Orientation]'];
        end
        obj                 = [obj; cell2mat(obj_tmp)]; %#ok<AGROW>
    else
        for jj = 1:num_overlap
            % slice binary radargram
            data_slice      = data_bin;
            data_slice(:, (jj * ind_step):ind_sep:end) ...
                            = false;
            % measure object properties of cleaned, sliced binary radargram
            curr_obj        = regionprops(labelmatrix(bwconncomp(bwmorph(data_slice, 'clean'), 4)), 'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
            % concatenate properties
            obj             = [obj; [cat(1, curr_obj.Centroid) [curr_obj(:).Area]' [curr_obj(:).MajorAxisLength]' [curr_obj(:).MinorAxisLength]' [curr_obj(:).Orientation]']]; %#ok<AGROW>
        end
    end
end