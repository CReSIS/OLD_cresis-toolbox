function P=Qchipr2(nu,lambda,x,epsilon)
%
%  This program computes the right-tail probability  
%  of a central or noncentral chi-squared PDF.
%
%  Input Parameters:
%
%    nu      = Degrees of freedom (1,2,3,etc.)
%    lambda  = Noncentrality parameter (must be positive),
%              set = 0 for central chi-squared PDF
%    x       = Real scalar value of random variable
%    epsilon = maximum allowable error (should be a small number
%              such as 1e-5) due to truncation of the infinite sum
%
%  Output Parameters:
%
%    P       = right-tail probability or the probability that the 
%              random variable exceeds the given value (1 - CDF)
%
%  Verification Test Case:
%
%    The inputs nu=1, lambda=2, x=0.5, epsilon=0.0001
%    should produce P=0.7772.
%    The inputs nu=5, lambda=6, x=10, epsilon=0.0001
%    should produce P=0.5063.
%    The inputs nu=8, lambda=10, x=15, epsilon=0.0001
%    should produce P=0.6161.
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

%  Determine how many terms in sum to be used (find M).
   t=exp(lambda/2)*(1-epsilon);
   sum=1;
   M=0;
   while sum < t
     M=M+1;
     sum=sum+((lambda/2)^M)/prod(1:M);
   end
%  Use different algorithms for nu even or odd.
   if (nu/2-floor(nu/2)) == 0  % nu is even.
%  Compute k=0 term of sum.
%  Compute Qchi2_nu(x).
%  Start recursion with Qchi2_2(x).
     Q2=exp(-x/2); g=Q2;
     for m=4:2:nu  %  If nu=2, loop will be omitted.
       g=g*x/(m-2); 
       Q2=Q2+g;
     end
%  Finish computation of k=0 term.
     P=exp(-lambda/2)*Q2;
%  Compute remaining terms of sum.
     for k=1:M
       m=nu+2*k;
       g=g*x/(m-2); Q2=Q2+g;
       arg=(exp(-lambda/2)*(lambda/2)^k)/prod(1:k);
       P=P+arg*Q2;
     end
   else  % nu is odd.
%  Compute k=0 term of sum.
     P=2*0.5*erfc(sqrt(x)/sqrt(2)); % Modified this line to avoid using "Q" function
%  Start recursion with Qchi2p_3(x).
     Q2p=sqrt(2*x/pi)*exp(-x/2); g=Q2p;
     if nu >1
       for m=5:2:nu  %  If nu=3, loop will be omitted.
         g=g*x/(m-2); 
         Q2p=Q2p+g;
       end
       P=P+exp(-lambda/2)*Q2p;
%  Compute remaining terms of sum.
       for k=1:M
         m=nu+2*k;
         g=g*x/(m-2); Q2p=Q2p+g;
         arg=(exp(-lambda/2)*(lambda/2)^k)/prod(1:k);
         P=P+arg*Q2p;
       end
     else
%  If nu=1, the k=0 term is just Qchi2_1(x)=2Q(sqrt(x)).
%  Add the k=0 and k=1 terms.
       P=P+exp(-lambda/2)*(lambda/2)*Q2p;
%  Compute remaining terms.
       for k=2:M
         m=nu+2*k;
         g=g*x/(m-2); Q2p=Q2p+g;
         arg=(exp(-lambda/2)*(lambda/2)^k)/prod(1:k);
         P=P+arg*Q2p;
       end
     end
   end    

