function I = Simpson2( f,a,b,N )
%
% Simpson estimates the value of the integral of f(x)
% from a to b by using the composite Simpson¡¯s 3/8 rule
% applied to n equal-length subintervals.
%
% I = Simpson(f,a,b,n) where
%
% f is an inline function representing the integrand,
% a, b are the limits of integration,
% 3N is the number of subintervals_,
%
% I is the integral estimate.
n=3*N;
h=(b-a)/n;
x=a:h:b;
I = 3*h/8*(2*sum(f(x(1:3:end))) + ...
    3*sum(f(x(2:3:end)))+ ...
    3*sum(f(x(3:3:end)))-f(a)-f(b));
end

