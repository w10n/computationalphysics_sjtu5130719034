function [r cc] = Secant2(f, x1, x2, tol, N)
%
% Secant uses secant method to approximate roots of
% f(x) = 0.
%
% [r n] = Secant(f, x1, x2, tol, N), where
%
% f is an inline function which represents f(x),
% x1 and x2 are the initial values of x,
% tol is the scalar tolerance of convergence (default
% is 1e-4),
% N is the maximum number of iterations (default is 20),
%
% r is the approximate root of f(x) = 0,
% n is the number of iterations required for
% convergence.
%
% Ramin S. Esfandiari, Numerical methods for engineers and scientists using Matlab
% Section 3.5, p. 83
cc=zeros(1,50);
if nargin < 5, N = 20;end
if nargin < 4,tol = 1e-4;end
x = zeros(1, N+1); % Pre-allocate
for n = 2:N,
    if x1 == x2
        r = 'failure';
        return
    end
    x(1) = x1; x(2) = x2;
    x(n+1) = x(n)-((x(n)-x(n-1))/(f(x(n))-f(x(n-1))))...
        * f(x(n));
    cc(n-1)=x(n+1);
    if abs(x(n+1)-x(n)) < tol
        r = x(n+1);
        return
    end
end
r = 'failure';