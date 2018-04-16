function I = TrapComp(f,a,b,n)
%
% TrapComp estimates the value of the integral of f(x)
% from a to b by using the composite trapezoidal rule
% applied to n equal-length subintervals.
%
% I = TrapComp(f,a,b,n) where
%
% f is an inline function representing the integrand,
% a and b are the limits of integration,
% n is the number of equal-length subintervals in [a,b],
%
% I is the integral estimate.
%
% Ramin S. Esfandiari, Numerical Methods for Engineers and Scientists Using
% Matlab,
% Section 6.3.4, p.291
%
h=(b-a)/n;
x=a:h:b;
I = h*(f(a)/2.+sum(f(x(2:n)))+f(b)/2);