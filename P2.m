close;
clear;
clc;

E=[0,25,50,75,100,125,150,175,200];
fE=[10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7];
e=[9.34,17.9,41.5,85.5,51.5,21.5,10.8,6.29,4.14];

fd=@(a) [f1(e,a,E,fE);f2(e,a,E,fE);f3(e,a,E,fE)];
r = M_N([60500 78 55^2/4]',fd,3,1e-4);
fr=r(1);ER=r(2);T=sqrt(r(3)*4);
display(fr);
display(ER);
display(T);

figure(1)
plot(E,fE,'ro');
hold on
x=0:1:200;
fEp=fr./((x-ER).^2+r(3));
plot(x,fEp);
hold off
legend('Experimental values','Fiting curve'),xlabel('E(MeV)'),ylabel('f(E)(mb)'),title('BW formula fit')
fprintf('Q1(d):Using BW formula,(fr,Er,T)=(%2.3d,%2.3d,%2.3d)\n',fr,ER,T)

function [xv,it] = M_N(x,f,n,tol)
% Newton Secant method for solving a system of n nonlinear equations
% in n variables.
% Example call: [xv,it] = Secantmv(x,f,n,tol)
% Requires an initial approximation column vector x. tol is
% required accuracy. User must define functions f (system equations)
% as a column vector. xv is the solution vector, the it
% parameter is number of iterations taken.
% WARNING. The method may fail, for example if initial estimates are poor.
%
it = 0; xv = x;
fr = feval(f,xv);
 while norm(fr) > tol
    Jr = jacobian(xv,f); xv = xv-Jr\fr;
    fr = feval(f,xv); it = it+1;
 end
 return
end

function [ output ] = f0( e,a,x )
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
output=a(1)/((x-a(2))^2+a(3));

end

function [ output ] = f1(e,a,x,y )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
output=0;
for n=1:9
    output=output+(y(n)-f0(e,a,x(n)))/(e(n)^2*((x(n)-a(2))^2+a(3)));
end
    
end

function [ output ] = f2(e,a,x,y )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
output=0;
for n=1:9
    output=output+(y(n)-f0(e,a,x(n)))*(x(n)-a(2))/(e(n)^2*((x(n)-a(2))^2+a(3)));
end

end

function [ output ] = f3(e,a,x,y )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
output=0;
for n=1:9
    output=output+(y(n)-f0(e,a,x(n)))/(e(n)^2*((x(n)-a(2))^2+a(3))^2);
end
    
end

function [J]=jacobian(x,f)

N=length(x);
J=zeros(N,N);
deltax=0.1*x;
newx=repmat(x',N,1);
for n=1:N
    newx(n,n)=newx(n,n)+deltax(n);
end

fv1=feval(f,x);
 
for n=1:N
    fv2=feval(f,newx(n,:));
    J(:,n)=(fv2-fv1)/deltax(n);
end
end