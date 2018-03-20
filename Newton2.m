function [r cc] = Newton2(f, fp, x1, tol, N)
%
% Newton uses Newton¡¯s method to approximate a root of
% f(x) = 0.
%
% [r n] = Newton(f, fp, x1, tol, N), where
%
% f is an inline function representing f(x),
% fp is an inline function representing f'(x),
% x1 is the initial point,
% tol is the scalar tolerance for convergence
% (default 1e-4),
% N is the maximum number of iterations (default 20),
%
% r is the approximate root of f(x) = 0,
% n is the number of iterations required for convergence.
% Ramin S. Esfandiari, Numerical methods for engineers and scientists using Matlab
% Section 3.5, p. 75
cc=zeros(1,50);
if nargin < 5, N = 20;end
if nargin < 4,tol = 1e-4;end
x = zeros(1, N+1); % Pre-allocate
x(1) = x1;
for n = 1:N,
    if fp(x(n)) == 0
        r = 'failure';
    return
    end
    x(n+1) = x(n)-f(x(n))/fp(x(n));
    cc(n)=x(n+1);
    if abs(x(n+1)-x(n)) < tol,
        r = x(n+1);
        return
    end
end
r = 'failure';

