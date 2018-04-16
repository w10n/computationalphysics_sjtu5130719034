function I = Simpson(f,a,b,n)
%
% Simpson estimates the value of the integral of f(x)
% from a to b by using the composite Simpson¡¯s 1/3 rule
% applied to n equal-length subintervals.
%
% I = Simpson(f,a,b,n) where
%
% f is an inline function representing the integrand,
% a, b are the limits of integration,
% n is the (even) number of subintervals,
%
% I is the integral estimate.
%
% Ramin S. Esfandiari, Numerical Methods for Engineers and Scientists Using
% Matlab,
% Section 6.3.5, p.295
%
h=(b-a)/n;
x=a:h:b;
I = h/3*(2*sum(f(x(1:2:end))) + ...
    4*sum(f(x(2:2:end)))-f(a)-f(b));