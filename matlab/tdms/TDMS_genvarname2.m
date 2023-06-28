function varName = TDMS_genvarname2(string,replaceStr,prependStr,alwaysPrepend)
%GENVARNAME2: Outputs a valid variable name VARNAME
%
%   GENVARNAME2(STRING,*REPLACESTR,*PREPENDSTR,*ALWAYSPREPEND) replaces invalid
%       characters in STRING with REPLACESTR. If the first letter is not
%       a letter, or ALWAYSPREPEND is true, then PREPENDSTR is prepended
%       to the variable name
%
%   INPUTS
%   =======================================================================
%   STRING         - input string to convert to safe variable name
%   *REPLACESTR    - (default '_'), character to replace invalid values with
%   *PREPENDSTR    - (default 'v'), value to prepend if first character is not
%                    a letter or ALWAYSPREPEND is true
%   *ALWAYSPREPEND - (default true), if true, alwyays adds the PREPENDSTR
%
%   EXAMPLES
%   ========================================
%   Example 1:
%       % always append is selected
%       % 'var' gets appended
%       % underscore is 
%       varName = genvarname2('RZ(1)','_','var',1)
%       varName = varRZ_1_
%   Example 2:
%       %always append is not selected
%       %first character is not numeric
%       varName = genvarname2('RZ(1)','_','var',0)
%       varName = RZ_1_
%   Example 3:
%       %always appedn is not selected
%       %since numeric appends 'var'
%       varName = genvarname2('1rack','_','var',0)
%       varName = var1rack
%
%   See also: genvarname
%
%   Copied from local genvarname2

if nargin < 2, replaceStr = []; end
if nargin < 3, prependStr = []; end
if nargin < 4, alwaysPrepend = []; end

if isempty(replaceStr), replaceStr = '_'; end
if isempty(prependStr), prependStr = 'v'; end
if isempty(alwaysPrepend), alwaysPrepend = true; end

if ~isempty(find(~(isstrprop(replaceStr,'alphanum') | replaceStr == '_'),1))
    error('REPLACESTR must be alphaNumeric or an underscore')
end

mask = isstrprop(string,'alphanum');
ind  = find(mask,1,'first');
string(~mask) = replaceStr; %Replaces all non-alpha numeric values with _
string = string(ind:end);

if(~isletter(string(1))) || alwaysPrepend
    if ~isempty(find(~(isstrprop(prependStr,'alphanum') | replaceStr == '_'),1))
        error('PREPENDSTR values must be alphaNumeric or an underscore')
    elseif ~isletter(prependStr(1))
        error('PREPENDSTR(1) needs to be a letter to have a valid variable name')
    end
    varName = [prependStr string];
else
    varName = string;
end