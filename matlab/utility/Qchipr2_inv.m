function x=Qchipr2_inv(nu,lambda,P,epsilon)
%
%  This program computes the inverse of Qchipr2.
%
%  Input Parameters:
%
%    nu      = Degrees of freedom (1,2,3,etc.)
%    lambda  = Noncentrality parameter (must be positive),
%              set = 0 for central chi-squared PDF
%    P       = right-tail probability or the probability that the 
%              random variable exceeds the given value (1 - CDF)
%    epsilon = maximum allowable error (should be a small number
%              such as 1e-5) due to truncation of the infinite sum
%
%  Output Parameters:
%
%    x       = Real scalar value of random variable
%
%  Verification Test Case:
%
%    The inputs nu=1, lambda=2, P=0.7772, epsilon=0.00001
%    should produce x=0.5.
%      x=Qchipr2_inv(1,2,0.7772,0.00001)
%    The inputs nu=5, lambda=6, P=0.5063, epsilon=0.00001
%    should produce x~=10.
%      x=Qchipr2_inv(5,6,0.5063,0.00001)
%    The inputs nu=8, lambda=10, P=0.6161, epsilon=0.00001
%    should produce x~=15.
%      x=Qchipr2_inv(8,10,0.6161,0.00001)
%
%
% From Mathworks Exchange:
% This directory contains the M-files described in the text 
% "Fundamentals of Statistical Signal Processing: Detection 
% Theory", by Steven Kay, Prentice Hall, 1998.  The files are:
% 
% dp.m
% montecarlo.m
% plotprob.m
% Q.m
% Qinv.m
% Qchipr2.m
% 
% They are provided free to the user and may be used without 
% restriction.   Their use in the preparation of any publication, 
% however, should be referenced.  The files were created using 
% MATLAB version 4.2 but should run under MATLAB versions 5.1, 5.2.  
% No toolboxes are required. 
% 
% Please note that:  These M-files are User Contributed Routines 
% which are being redistributed by the MathWorks, upon request, on 
% an "as is" basis.  A User Contributed Routine is not a product of 
% The Mathworks, Inc. and The Mathworks, Inc. assumes no 
% responsibility for any errors that may exist in these routines.
% 
% Any questions on this software should be directed to:
% 
% Professor Steven Kay
% Dept. of Electrical and Computer Engineering
% Kelley Hall
% 4 East Alumni Ave.
% University of Rhode Island
% Kingston, RI 02881
% 
% phone: 401-874-5804
% fax: 401-782-6422
% email:kay@ele.uri.edu

x = fminsearch(@(x) abs(Qchipr2(nu,lambda,x,epsilon) - P),1);
