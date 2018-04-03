function c = Bisection(f, a, b, kmax, tol)
%
% Bisection uses the bisection method to approximate a
% root of f(x) = 0 in the interval [a,b].
%
% c = Bisection(f, a, b, kmax, tol) where
%
% f is an inline function representing f(x),
% a and b are the limits of the interval [a,b],
% kmax is the maximum number of iterations
% (default 20),
% tol is the scalar tolerance for convergence
% (default 1e¨C4),
%
% c is the approximate root of f(x) = 0.
%
% Ramin S. Esfandiari, Numerical methods for engineers and scientists using Matlab
% Section 3.2, p. 60
cc=zeros(1,50);
if nargin < 5, tol = 1e-4; end
if nargin < 4, kmax = 20; end
if f(a)*f(b) > 0
    c ='failure';
    return
end
%disp(' k a b c (b¨Ca)/2')
for k = 1:kmax
    c =(a+b)/2; % Find the first midpoint
    cc(k)=c;
    if f(c) == 0 % Stop if a root has been found
        return
    end
    %fprintf('%3i %11.6f%11.6f%11.6f%11.6f\n',k,a,b,c,(b-a)/2)
    if (b-a)/2 < tol % Stop if tolerance is met
        return
    end
    if f(b)*f(c) > 0 % Check sign changes
        b = c; % Adjust the endpoint of interval
    else
        a = c;
    end
end
c = 'failure';
