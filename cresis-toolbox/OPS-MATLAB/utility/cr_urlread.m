function [output,status] = opsUrlRead(urlChar,user, pass, method, params)
% CReSIS urlread function supporting basic authentication. 
% 
% Input: 
%   urlChar: string of url
%   user: string - username for authentication 
%   pass: string - password for authentication 
%   method: string - 'POST' or 'GET'
%   params: structure - {'field_nm' data}
%
% Output: 
%   output: string - html elements of url and any return
%   status: optional. 
%
% See also: urlread.m, urlreadwrite.m
% Authors: Trey Stafford and Kyle Purdon

if ~usejava('jvm')
   error(message('MATLAB:urlread:NoJvm'));
end

import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;

% Be sure the proxy settings are set.
com.mathworks.mlwidgets.html.HTMLPrefs.setProxySettings

% Check number of inputs and outputs.
error(nargchk(3,5,nargin))
error(nargoutchk(0,2,nargout))
if ~ischar(urlChar)
    error('MATLAB:urlread:InvalidInput','The first input, the URL, must be a character array.');
end
if (nargin > 3) && ~strcmpi(method,'get') && ~strcmpi(method,'post')
    error('MATLAB:urlread:InvalidInput','Second argument must be either "get" or "post".');
end

% Do we want to throw errors or catch them?
if nargout == 2
    catchErrors = true;
else
    catchErrors = false;
end

% Set default outputs.
output = '';
status = 0;

% Encode username and password to base64
passstr = horzcat(user, ':', pass);
encoder = sun.misc.BASE64Encoder();
out = char(encoder.encode(java.lang.String(passstr).getBytes()));
auth = horzcat('Basic ', out);


% GET method.  Tack param/value to end of URL.
if (nargin > 3) && strcmpi(method,'get')
    if mod(length(params),2) == 1
        error('MATLAB:urlread:InvalidInput','Invalid parameter/value pair arguments.');
    end
    for i=1:2:length(params)
        if (i == 1), separator = '?'; else separator = '&'; end
        param = char(java.net.URLEncoder.encode(params{i}));
        value = char(java.net.URLEncoder.encode(params{i+1}));
        urlChar = [urlChar separator param '=' value];
    end
end

% Create a urlConnection.
[urlConnection,errorid,errormsg] = urlreadwrite(mfilename,urlChar);
if isempty(urlConnection)
    if catchErrors, return
    else error(errorid,errormsg);
    end
end

%Use encoded username/password to authenticate
urlConnection.setRequestProperty(...
    'Authorization', auth);

% POST method.  Write param/values to server.
if (nargin > 3) && strcmpi(method,'post')
    try
        urlConnection.setDoOutput(true);
        urlConnection.setRequestProperty( ...
            'Content-Type','application/x-www-form-urlencoded');
        printStream = java.io.PrintStream(urlConnection.getOutputStream);
        for i=1:2:length(params)
            if (i > 1), printStream.print('&'); end
            param = char(java.net.URLEncoder.encode(params{i}));
            value = char(java.net.URLEncoder.encode(params{i+1}));
            printStream.print([param '=' value]);
        end
        printStream.close;
    catch
        if catchErrors, return
        else error('MATLAB:urlread:ConnectionFailed','Could not POST to URL.');
        end
    end
end

% Read the data from the connection.
try
    inputStream = urlConnection.getInputStream;
    byteArrayOutputStream = java.io.ByteArrayOutputStream;
    % This StreamCopier is unsupported and may change at any time.
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    isc.copyStream(inputStream,byteArrayOutputStream);
    inputStream.close;
    byteArrayOutputStream.close;
    output = char(byteArrayOutputStream.toByteArray');
catch
    if catchErrors, return
    else error('MATLAB:urlread:ConnectionFailed','Error downloading URL. Your network connection may be down or your proxy settings improperly configured.');
    end
end

status = 1;