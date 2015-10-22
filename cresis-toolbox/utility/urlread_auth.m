function [s,info] = urlread_auth(url, user, password)
%URLREAD_AUTH Like URLREAD, with basic authentication
%
% [s,info] = urlread_auth(url, user, password)
%
% Returns bytes. Convert to char if you're retrieving text.
%
% Examples:
% sampleUrl = 'http://browserspy.dk/password-ok.php';
% [s,info] = urlread_auth(sampleUrl, 'test', 'test');
% txt = char(s)
%

% Matlab's urlread() doesn't do HTTP Request params, so work directly with Java
jUrl = java.net.URL(url);
conn = jUrl.openConnection();
conn.setRequestProperty('Authorization', ['Basic ' base64encode([user ':' password])]);
conn.connect();
info.status = conn.getResponseCode();
info.errMsg = char(readstream(conn.getErrorStream()));
s = readstream(conn.getInputStream());

function out = base64encode(str)
% Uses Sun-specific class, but we know that is the JVM Matlab ships with
encoder = sun.misc.BASE64Encoder();
out = char(encoder.encode(java.lang.String(str).getBytes()));

%%
function out = readstream(inStream)
%READSTREAM Read all bytes from stream to uint8
try
    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;
    byteStream = java.io.ByteArrayOutputStream();
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier();
    isc.copyStream(inStream, byteStream);
    inStream.close();
    byteStream.close();
    out = typecast(byteStream.toByteArray', 'uint8'); %'
catch err
    out = []; %HACK: quash
end