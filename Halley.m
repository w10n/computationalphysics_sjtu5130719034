function [c cc] = Halley(f,fp,fpp, x1, kmax, tol)
%using Halley's method to find to roots of certain function
% root of f(x) = 0 near a
% fp is the derivative of f
% fpp is the derivative of fp
%c is the root founded
%kmax is the maximum number of iterations
%(default 20),
% tol is the scalar tolerance for convergence
% (default 1e¨C4)
%cc is all the x calculated

cc=zeros(1,50);
if nargin < 5, tol = 1e-4; end
if nargin < 4, kmax = 20; end
cc=zeros(1,50);
x(1)=x1;
for n=1:kmax,
    x(n+1) =x(n)-2*f(x(n))*fp(x(n))/(2*fp(x(n))^2-f(x(n))*fpp(x(n)));
    cc(n)=x(n);
    if abs(x(n+1)-x(n)) < tol,
        c = x(n+1);
        return
    end
end
c = 'failure';
    
